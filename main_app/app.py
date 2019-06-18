import json
import os
import sys
import uuid
from pathlib import Path
from shutil import make_archive, rmtree

from flask import (Flask, flash, jsonify, redirect, render_template, request,
                   session, url_for)
from livereload import Server
from werkzeug.utils import secure_filename
from wtforms import (DecimalField, FieldList, FileField, FormField,
                     IntegerField, StringField)
from wtforms.validators import InputRequired, NumberRange

from static.python.analyze import (check_filetypes, check_num_conditions,
                                   create_experiment_record,
                                   generate_analyze_kwargs,
                                   get_experiment_name, run_spplot_analysis,
                                   upload_datafiles_to_tmp,
                                   write_config_to_tmp_data_upload_dir)
from static.python.plot import generate_plots
from static.python.utils import make_temp_upload_dir

os.environ['FLASK_DEBUG'] = '1'
os.environ['FLASK_APP'] = 'app.py'
os.environ['FLASK_ENV'] = 'development'

app = Flask(__name__)

MAX_SAMPLE_ROWS = 11
MIN_SAMPLE_ROWS = 1
SPPLOT_BUCKET = 'spplot'

SECRET_KEY = os.urandom(32)
UPLOAD_FOLDER = 'tmp'


app.config['SECRET_KEY'] = SECRET_KEY
ALLOWED_EXTENSIONS = ['.spp']
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


@app.route('/index')
@app.route('/')
def index():
    return render_template('index.html')


@app.route('/analyze', methods=['GET', 'POST'])
def analyze():

    if request.method == "POST":
        " if form is submitted using POST "
        unique_data_id = str(uuid.uuid4())
        tmp_data_upload_dir = os.path.join(app.config['UPLOAD_FOLDER'], unique_data_id)
        tmp_output_dir = os.path.join(tmp_data_upload_dir, 'analysis_results')
        os.makedirs(tmp_data_upload_dir)
        os.makedirs(tmp_output_dir)

        session['tmp_data_upload_dir'] = tmp_data_upload_dir
        session['tmp_output_dir'] = tmp_output_dir

        num_conditions = check_num_conditions(request)

        control_samples, treated_samples = upload_datafiles_to_tmp(request, num_conditions, tmp_data_upload_dir)
        if not control_samples:
            rmtree(tmp_data_upload_dir)
            return render_template('analyzing/analyze.html', status="File upload Failed")

        files_ok, fail_list = check_filetypes(tmp_data_upload_dir, control_samples, treated_samples)
        if not files_ok:
            rmtree(tmp_data_upload_dir)
            return render_template('analyzing/analyze.html', status="File file type check failed on: {}".format(fail_list))

        spplot_filter_kwargs = generate_analyze_kwargs(
            request,
            num_conditions,
            control_samples,  # files names
            treated_samples,  # file file_names
            tmp_data_upload_dir,
            tmp_output_dir)


        result, err = run_spplot_analysis(**spplot_filter_kwargs)


        if result is False:
            print("RESULT FAILED")
            return render_template('analyzing/analyze.html', status='{}'.format(err))

        # for rerunning locally
        config_string = create_experiment_record(**spplot_filter_kwargs)
        config_path = write_config_to_tmp_data_upload_dir(config_string, tmp_data_upload_dir)

        new_dir = os.path.join('static', 'analysis_results', unique_data_id)
        os.makedirs(new_dir)
        print(new_dir)
        # make archive to send to s3 for download
        output_path = os.path.join(app.config['UPLOAD_FOLDER'], "_".join([unique_data_id[:13], get_experiment_name(request)]))
        make_archive(
            output_path,
            'zip',
            app.config['UPLOAD_FOLDER']
            )
        os.rename(tmp_data_upload_dir, new_dir)

        print("UPLOADING TO S3 FAKE")
        # upload_to_s3(output_path)  # need to make sure the right creds are used to do this with my account
        os.remove(output_path + '.zip')
        image_files = [os.path.join(new_dir, 'analysis_results', x) for x in os.listdir(os.path.join(new_dir, 'analysis_results')) if x.endswith('png')]

        return render_template('analyzing/analysis_result.html', status="SUCCESS!", image_files=enumerate(image_files), link='https://www.google.com')

    return render_template('analyzing/analyze.html', status='Experiment not yet run.')#, param_form=param_form)


@app.route('/interpretting')
def interpretting():
    return render_template('interpretting.html')


@app.route('/make_plots', methods=['GET', 'POST'])
def plotting():

    if request.method == 'POST':

        unique_data_id = str(uuid.uuid4())
        tmp_data_upload_dir = session.get('tmp_data_upload_dir') or os.path.join(app.config['UPLOAD_FOLDER'], unique_data_id)
        tmp_output_dir = session.get('tmp_output_dir') or os.path.join(tmp_data_upload_dir, 'plotting_results')

        os.makedirs(tmp_data_upload_dir)
        os.makedirs(tmp_output_dir)

        data_type = 'plot'
        image_files = generate_plots(
            request,
            data_type,
            tmp_output_dir,
            min_diff=0,
            max_diff=10000
        )

        new_dir = os.path.join('static', 'plot_results', unique_data_id)
        os.makedirs(new_dir)
        static_image_files = list()
        for old_path in image_files:
            new_path = os.path.join(new_dir, os.path.basename(old_path))
            os.rename(old_path, new_path)
            static_image_files.append(new_path)

        return render_template('plotting/make_plots.html', image_files=enumerate(static_image_files))

    return render_template('plotting/make_plots.html', image_files=enumerate(['']))


@app.route('/troubleshooting')
def troubleshooting():
    return render_template('troubleshooting.html')


@app.route('/citing')
def citing():
    return render_template('citing.html')


if __name__ == '__main__':
    make_temp_upload_dir(UPLOAD_FOLDER)
    app.run()
