import json
import os
import sys
import uuid
from shutil import rmtree, make_archive

from flask import (Flask, flash, jsonify, redirect, render_template, request,
                   url_for, session)
from flask_wtf import FlaskForm
from jinja2 import Template
from livereload import Server
from werkzeug.utils import secure_filename
from wtforms import (DecimalField, FieldList, FileField, FormField,
                     IntegerField, StringField)
from wtforms.validators import InputRequired, NumberRange

import boto3
from main_app.seqpyplot.analyzer.paired_sample_filter import PairedSampleFilter
from main_app.seqpyplot.container.data_container import DataContainer
from main_app.seqpyplot.parsers.config_parser import config_parser
from main_app.seqpyplot.parsers.gene_list_parser import MakeFigureList
from main_app.seqpyplot.plot.bar_plotter import PairedBarPlot
from main_app.seqpyplot.plot.de_tally_plotter import TallyDe
from main_app.seqpyplot.plot.paired_line_plotter import PairedDataLinePlotter
from main_app.seqpyplot.plot.PCA import PCADecomposition
from main_app.seqpyplot.plot.scatter_plotter import ScatterPlots
from main_app.seqpyplot.printers.data_printer import DataPrinter
from main_app.seqpyplot.utils import make_default_output_dir

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
ALLOWED_EXTENSIONS = ['.splot']
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


if os.path.exists(UPLOAD_FOLDER):
    rmtree(UPLOAD_FOLDER)
    os.mkdir(UPLOAD_FOLDER)
else:
    os.mkdir(UPLOAD_FOLDER)


def check_files_present(request_files):
    if len(list(request_files.keys())) == 0:
        return False
    return True


def check_suffices(request_files):
    if len([x for x in request_files.keys()]) % 2 == 0:
        if all([x.filename.endswith('.splot') for x in request_files.values()]):
            return True
    return False


def check_filetypes(tmp_dir, control_samples, treated_samples):
    "check all files first line is elements and last one can be int"
    files_check = True
    fail_list = list()
    files = control_samples + treated_samples
    for file_ in files:
        with open(os.path.join(tmp_dir, file_), 'r') as fin:
            first_line = fin.readline()
            split = first_line.split()

            try:
                int(split[-1]) # last el should be int
            except ValueError:
                fail_list.append(file_)
                files_check = False

            if not len(split) == 2:
                fail_list.append(file_)
                files_check = False

    return files_check, fail_list


def check_num_conditions(request):
    conditions = set([x.split('_')[0] for x in request.files])
    return len(conditions)


def upload_datafiles_to_tmp(request, num_conditions, tmp_dir):

    files = [x for x in request.files]

    control_samples = list()
    treated_samples = list()

    if num_conditions == 1:

        for file_ in request.files:
            name = request.files[file_].name
            filename = request.files[file_].filename
            request.files[file_].save(os.path.join(tmp_dir, filename))

            control_samples.append(filename)

    elif num_conditions == 2:
        for file_ in request.files:
            name = request.files[file_].name
            filename = request.files[file_].filename
            request.files[file_].save(os.path.join(tmp_dir, filename))

            if 'control' in file_:
                control_samples.append(filename)
            elif 'treated' in file_:
                treated_samples.append(filename)
    else:
        return False, False

    return control_samples, treated_samples


def retrieve_form_values(form, value_type):
    keys = [x for x in filter(lambda x: value_type in x, form)]
    values = [form.get(key) for key in keys]
    return keys, values


def get_experiment_name(request):
    experiment_ids, experiment_name = retrieve_form_values(form=request.form, value_type='experimentname')
    assert len(experiment_name) == 1
    return experiment_name.pop()


def get_time_point_names(request):
    time_point_ids, time_point_names = retrieve_form_values(form=request.form, value_type='time')
    return time_point_names


def get_param_dict(request):

    param_ids, param_values = retrieve_form_values(form=request.form, value_type='param')
    param_names = [  # names need to match the analysis config template jinja
        'log2fold', 'expression_min', 'expression_max', 'min_diff', 'max_diff',
        'correction_interval_min', 'correction_interval_max', 'num_components'
        ]

    from collections import OrderedDict
    param_dict = dict()
    for name, value in zip(param_names, param_values):
        param_dict[name] = safe_string_to_number(value)
    return param_dict


def safe_string_to_number(number):
    if '.' in number:
        try:
            number = float(number)
        except ValueError:
            pass
    else:
        try:
            number = int(number)
        except ValueError:
            pass
    return number


def render_config_template(
    counts_dir,
    analysis_output_dir,
    controls_list,
    treated_list,
    data_type,
    control_sample_names,
    treated_sample_names,
    time_point_names,
    condition_names,
    experiment_name,
    log2fold,
    expression_min,
    expression_max,
    min_diff,
    max_diff,
    correction_interval_min,
    correction_interval_max,
    num_components):

    config_template = Template(
"""[data_directory]
dir={{ counts_dir }}
output={{ analysis_output_dir }}

[data]
controls={{ controls_list }}
treated={{ treated_list }}
data_type={{ data_type }}

[names]
controls={{ control_sample_names }}
treated={{ treated_sample_names }}
times={{ time_point_names }}
conditions={{ condition_names }}
experiment_name={{ experiment_name }}

[params]
log2fold={{ log2fold }}
low={{ expression_min }}
hi={{ expression_max }}
diff=[{{ min_diff }}, {{ max_diff }}]
correction_interval_min={{ correction_interval_min }}
correction_interval_max={{ correction_interval_max }}
num_components={{ num_components }}

[file_names]
prefix={{ experiment_name }}

[plot_options]
scatrange=[10, 10000]
""")

    analysis_config = config_template.render(
        counts_dir = counts_dir,
        analysis_output_dir = analysis_output_dir,
        controls_list = controls_list,
        treated_list = treated_list,
        data_type = data_type,
        control_sample_names= control_sample_names,
        treated_sample_names = treated_sample_names,
        time_point_names = time_point_names,
        condition_names = condition_names,
        experiment_name = experiment_name,
        log2fold = log2fold,
        expression_min = expression_min,
        expression_max = expression_max,
        min_diff = min_diff,
        max_diff = max_diff,
        correction_interval_min = correction_interval_min,
        correction_interval_max = correction_interval_max,
        num_components = num_components)

    return analysis_config


def generate_template_kwargs(request, num_conditions, control_samples, treated_samples, tmp_dir, tmp_output_dir):
    template_kwargs = dict()

    template_kwargs['counts_dir'] = tmp_dir
    template_kwargs['analysis_output_dir'] = tmp_output_dir

    template_kwargs['controls_list'] = control_samples
    template_kwargs['treated_list'] = treated_samples

    data_type = get_datatype(request)
    template_kwargs['data_type'] = data_type

    template_kwargs['control_sample_names'] = [x.split('.')[0] for x in control_samples]
    template_kwargs['treated_sample_names'] = [x.split('.')[0] for x in treated_samples]

    time_point_names = get_time_point_names(request)
    template_kwargs['time_point_names'] = time_point_names

    if num_conditions == 1:
        template_kwargs['condition_names'] = ['single_condition']
    elif num_conditions == 2:
        template_kwargs['condition_names'] = ['control', 'treated']
    else:
        template_kwargs['condition_names'] = [str(x) for x in range(len(num_conditions))]

    experiment_name = get_experiment_name(request)  # string
    template_kwargs['experiment_name'] = experiment_name

    params_dict =  get_param_dict(request)  # dictionary
    template_kwargs.update(params_dict)

    return template_kwargs


def write_config_to_tmp_dir(config_string, tmp_dir):
    config_path = os.path.join(tmp_dir, 'analysis_config')
    with open(config_path, 'w+') as fout:
        fout.write(config_string)

    return config_path


def get_datatype(request):
    # TODO collect from form?
    return 'htseq'


def upload_to_s3(file_to_upload):
    s3 = boto3.client('s3', region='us-east-2')
    key = os.path.join('SPLOT_USER_EXPERIMENTS', os.path.basename(file_to_upload))
    s3.upload_file(
        file_to_upload,
        SPPLOT_BUCKET,
        key,
        ExtraArgs={
            'ACL': 'public-read',
            'ContentType': "application/zip"
        }
    )

@app.route('/index')
@app.route('/')
def index():
    return render_template('index.html')


@app.route('/analyze', methods=['GET', 'POST'])
def analyze():

    if request.method == "POST":
        " if form is submitted using POST "
        unique_data_id = str(uuid.uuid4())
        tmp_dir = os.path.join(app.config['UPLOAD_FOLDER'], unique_data_id)
        tmp_output_dir = os.path.join(tmp_dir, 'analysis_results')
        os.makedirs(tmp_dir)
        os.makedirs(tmp_output_dir)

        num_conditions = check_num_conditions(request)

        control_samples, treated_samples = upload_datafiles_to_tmp(request, num_conditions, tmp_dir)
        if not control_samples:
            rmtree(tmp_dir)
            return render_template('analyzing/analyze.html', status="File upload Failed")

        files_ok, fail_list = check_filetypes(tmp_dir, control_samples, treated_samples)
        if not files_ok:
            rmtree(tmp_dir)
            return render_template('analyzing/analyze.html', status="File file type check failed on: {}".format(fail_list))

        template_kwargs = generate_template_kwargs(
            request,
            num_conditions,
            control_samples,
            treated_samples,
            tmp_dir,
            tmp_output_dir)

        config_string = render_config_template(**template_kwargs)
        config_path = write_config_to_tmp_dir(config_string, tmp_dir)

        result, err = run_spplot_analysis(config_path)
        if result is False:
            print("RESULT FAILED")
            return render_template('analyzing/analyze.html', status='{}'.format(err))

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
        os.rename(tmp_dir, new_dir)

        print("UPLOADING TO S3 FAKE")
        # upload_to_s3(output_path)  # need to make sure the right creds are used to do this with my account
        os.remove(output_path + '.zip')

        session['session_data'] = {
            'result_images': [x for x in os.listdir(os.path.join(new_dir, 'analysis_results')) if x.endswith('png')],
            'unique_key': unique_data_id
        }
        return redirect(url_for('success'))

    return render_template('analyzing/analyze.html', status='Experiment not yet run.')#, param_form=param_form)


@app.route('/success', methods=['GET'])  #TODO REMOVE GET as method -- this should only be available AFTER analysis so by redirect only
def success():
    session_key = session.get('session_data', {}).get('unique_key', None)
    files = session.get('session_data', {}).get('result_images', None)
    image_files = list(enumerate([os.path.join('static', 'analysis_results', session_key, 'analysis_results', x) for x in files]))
    return render_template('analyzing/analysis_result.html', status="SUCCESS!", image_files=image_files, link='link')


@app.route('/interpretting')
def interpretting():
    return render_template('interpretting.html')


def check_num_series(request):
    number_map = {
        'one': 1,
        'two': 2,
        'three': 3,
        'four': 4,
        'five': 5
    }
    return number_map[request.values['num_series']]


def collect_gene_list(request):
    raw_list = request.values['gene_list']
    processed_list = [x.strip() for x in raw_list.split()]  # list of genes
    return processed_list


def collect_plot_config_path(request, tmp_plot_dir):
    name = request.files['analysis_config'].name
    filename = request.files['analysis_config'].filename

    file_path = os.path.join(tmp_plot_dir, filename)
    request.files['analysis_config'].save(file_path)
    return file_path


def load_de_genelist(file_path):
    with open(file_path, 'r') as fin:
        genes = [x.strip() for x in fin.readlines()]
    return genes


def upload_plottable_data(request, tmp_plot_dir):
    name = request.files['data_file'].name
    filename = request.files['data_file'].filename

    file_path = os.path.join(tmp_plot_dir, filename)
    request.files['data_file'].save(file_path)
    return file_path


def upload_de_genelist(request, tmp_plot_dir):
    name = request.files['de_gene_list'].name
    filename = request.files['de_gene_list'].filename

    file_path = os.path.join(tmp_plot_dir, filename)
    request.files['de_gene_list'].save(file_path)
    gene_list = load_de_genelist(file_path)
    return gene_list

def write_textbox_genelist_to_tmp(genelist, tmp_plot_dir):
    output_path = os.path.join(tmp_plot_dir, 'genelist.txt')
    with open(output_path, 'w') as fout:
        for gene in genelist:
            fout.write(gene + '\n')
    return output_path

@app.route('/make_plots', methods=['GET', 'POST'])
def plotting():

    if request.method == 'POST':
        unique_data_id = str(uuid.uuid4())
        tmp_dir = os.path.join(app.config['UPLOAD_FOLDER'], unique_data_id)
        tmp_output_dir = os.path.join(tmp_dir, 'plotting_results')
        os.makedirs(tmp_dir)
        os.makedirs(tmp_output_dir)

        num_conditions = check_num_series(request)
        processed_gene_list = collect_gene_list(request)
        gene_list_path = write_textbox_genelist_to_tmp(processed_gene_list, tmp_output_dir)
        config_path = collect_plot_config_path(request, tmp_output_dir)
        plottable_data_path = upload_plottable_data(request, tmp_output_dir)
        complete_de_genelist = upload_de_genelist(request, tmp_output_dir)

        config_obj = config_parser(config_path)
        config_obj.set('data_directory', 'dir', plottable_data_path)
        config_obj.set('data_directory', 'output', tmp_output_dir)
        config_obj.set('data', 'data_type', 'plot')
        config_obj.set('file_names', 'genelist', gene_list_path)

        # load the data container_obj
        container_obj = DataContainer(config_obj)
        data, ercc_data = container_obj.parse_plotter_data(plottable_data_path)
        filter_obj = PairedSampleFilter(config_obj)
        if complete_de_genelist:
            filter_obj.complete_de_gene_list = complete_de_genelist
        else:
            filter_obj.complete_de_gene_list = list()

        print("\nPlotting data...\n")
        line_plotter = PairedDataLinePlotter(config_obj, filter_obj, data)
        fig_list = MakeFigureList(config_obj)
        line_plotter.plot_figure(figure_list=fig_list.plot_groups, plottable_data=data)
        image_files = [os.path.join(tmp_output_dir, x) for x in os.listdir(tmp_output_dir) if x.endswith('png')]

        new_dir = os.path.join('static', 'plot_results', unique_data_id)
        os.makedirs(new_dir)
        static_image_files = list()
        for old_path in image_files:
            new_path = os.path.join(new_dir, os.path.basename(old_path))
            os.rename(old_path, new_path)
            static_image_files.append(new_path)

        # import pdb; pdb.set_trace()
        return render_template('plotting/make_plots.html', image_files=enumerate(static_image_files))


    return render_template('plotting/make_plots.html', image_files=enumerate(['']))


@app.route('/troubleshooting')
def troubleshooting():
    return render_template('troubleshooting.html')


@app.route('/citing')
def citing():
    return render_template('citing.html')


def run_spplot_analysis(config_path, correct_by_rotation=False):
    try:
        config_obj = config_parser(config_path)

        # load the data container_obj
        container_obj = DataContainer(config_obj)
        data, ercc_data = container_obj.parse_input()

        # if args.impute:
        #     print('Imputation not yet implemented')

        # TODO allow this option?
        # if not args.unnorm:
        data = container_obj.normalize_file_pairs(data) # Single df of normalized data

        split_data = container_obj.split(data)  # List of normalized dfs

        # TODO allow this?
        # if args.svd:
        #     split_data = container_obj.remove_variance(split_data)

        # TODO allow this?
        # if correct_by_rotation:
        #     split_data = container_obj.correct_via_rotation(split_data)

        #--------------------------------------------------------------------
        #  Filter data

        filter_obj = PairedSampleFilter(config_obj)
        filter_result = filter_obj.main_filter_process(split_data)

        #--------------------------------------------------------------------
        # Save filter results

        output_path = config_obj.get('data_directory', 'output')
        output_path = make_default_output_dir(output_path or None, overwrite=True)

        data_printer = DataPrinter(config_obj, container_obj or None, filter_obj or None)()

        #--------------------------------------------------------------------
        # Generate Plots

        # !!!TODO!!! seperate script here for plotting the data

        # print("\nPlotting data...\n")
        # line_plotter = PairedDataLinePlotter(config_obj, filter_obj, data)
        # fig_list = MakeFigureList(config_obj)
        # line_plotter.plot_figure(figure_list=fig_list.plot_groups, plottable_data=data)

        bar_plotter = PairedBarPlot(config_obj=config_obj)
        bar_plotter.create_bar_plot(filter_obj.de_count_by_stage)

        scatter_plotter = ScatterPlots(config_obj=config_obj, container_obj=container_obj, filter_obj=filter_obj)
        scatter_plotter.create_scatter_plots()

        tally_plotter = TallyDe(config_obj, container_obj)
        tally_plotter.create_tally_plot(split_data)

        pca_decomp = PCADecomposition(config_obj, container_obj)
        pca_decomp.create_pca_plot()

        #--------------------------------------------------------------------
        # Tidy up

        print("\nScript completed no errors")
        return True, ''

    except Exception as e:
        print("SCRIPT FKUP")
        return False, e


if __name__ == '__main__':
    if sys.argv[1] == 'debug':
        app.debug = True
        server = Server(app.wsgi_app)
        server.serve()
    else:
        app.run()
