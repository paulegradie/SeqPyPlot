from static.python.utils import retrieve_form_values
import os
from jinja2 import Template
from pathlib import Path


from seqpyplot.analyzer.paired_sample_filter import PairedSampleFilter
from seqpyplot.container.data_container import DataContainer
from seqpyplot.parsers.config_parser import config_parser
from seqpyplot.parsers.gene_list_parser import MakeFigureList
from seqpyplot.plot.bar_plotter import PairedBarPlot
from seqpyplot.plot.de_tally_plotter import TallyDe
from seqpyplot.plot.paired_line_plotter import PairedDataLinePlotter
from seqpyplot.plot.PCA import PCADecomposition
from seqpyplot.plot.scatter_plotter import ScatterPlots
from seqpyplot.printers.data_printer import DataPrinter

from static.python.utils import make_temp_upload_dir, safe_string_to_number


def get_experiment_name(request):
    experiment_ids, experiment_name = retrieve_form_values(form=request.form, value_type='experimentname')
    assert len(experiment_name) == 1
    return experiment_name.pop()


def get_datatype(request):
    # TODO collect from form?
    return 'htseq'


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


def get_time_point_names(request):
    time_point_ids, time_point_names = retrieve_form_values(form=request.form, value_type='time')
    return time_point_names


def generate_analyze_kwargs(request, num_conditions, control_samples, treated_samples, tmp_data_upload_dir, tmp_output_dir):
    template_kwargs = dict()

    template_kwargs['data_directory'] = tmp_data_upload_dir
    template_kwargs['output_dir'] = tmp_output_dir

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


def write_config_to_tmp_data_upload_dir(config_string, tmp_data_upload_dir):
    config_path = os.path.join(tmp_data_upload_dir, 'analysis_config')
    with open(config_path, 'w+') as fout:
        fout.write(config_string)

    return config_path


def create_experiment_record(
    data_directory,
    output_dir,
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
dir={{ data_directory }}
output={{ output_dir }}

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
expression_min={{ expression_min }}
expression_max={{ expression_max }}
min_diff={{ min_diff }}
max_diff={{ max_diff }}
correction_interval_min={{ correction_interval_min }}
correction_interval_max={{ correction_interval_max }}
num_components={{ num_components }}

[file_names]
prefix={{ experiment_name }}

[plot_options]
scatrange=[10, 10000]
""")

    analysis_config = config_template.render(
        data_directory = data_directory,
        output_dir = output_dir,
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
        num_components = num_components
        )

    return analysis_config


def check_filetypes(tmp_data_upload_dir, control_samples, treated_samples):
    "check all files first line is elements and last one can be int"
    files_check = True
    fail_list = list()
    files = control_samples + treated_samples
    for file_ in files:
        with open(os.path.join(tmp_data_upload_dir, file_), 'r') as fin:
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


def upload_datafiles_to_tmp(request, num_conditions, tmp_data_upload_dir):

    files = [x for x in request.files]

    control_samples = list()
    treated_samples = list()

    if num_conditions == 1:

        for file_ in request.files:
            name = request.files[file_].name
            filename = request.files[file_].filename
            request.files[file_].save(os.path.join(tmp_data_upload_dir, filename))

            control_samples.append(filename)

    elif num_conditions == 2:
        for file_ in request.files:
            name = request.files[file_].name
            filename = request.files[file_].filename
            request.files[file_].save(os.path.join(tmp_data_upload_dir, filename))

            if 'control' in file_:
                control_samples.append(filename)
            elif 'treated' in file_:
                treated_samples.append(filename)
    else:
        return False, False

    return control_samples, treated_samples


def run_spplot_analysis(
    data_directory,  # tmp data upload directory
    output_dir,  # tmp output directory
    controls_list,  # file name only
    treated_list,  # file name only
    data_type,
    control_sample_names,
    treated_sample_names,
    time_point_names,
    condition_names,
    experiment_name,

    # filter params
    log2fold,
    expression_max,
    expression_min,
    min_diff,
    max_diff,
    correction_interval_min,
    correction_interval_max,
    num_components,
    scatter_min=0,
    scatter_max=12000
):

    # prep the file paths
    data_directory = Path(data_directory)
    data_paths = [str(data_directory / x) for x in controls_list] + [str(data_directory / x) for x in treated_list] if treated_sample_names else []
    sample_names = control_sample_names + treated_sample_names
    file_name_pairs = [(x, y) for x, y in zip(control_sample_names, treated_sample_names)]
    num_file_pairs = len(file_name_pairs)
    # try:
    # load the data container_obj
    container_obj = DataContainer(
        data_directory=data_directory,
        data_paths=data_paths,
        sample_names=sample_names,
        time_point_names=time_point_names,
        file_name_pairs=file_name_pairs,
        num_file_pairs=num_file_pairs,
        num_components=num_components,
        data_type=data_type)

    # TODO alexpression_min for unnormed data to be submitted
    data, ercc_data = container_obj.parse_input()
    normalized_plottable_data = container_obj.normalize_file_pairs(data) # Single df of normalized data
    split_data = container_obj.split(normalized_plottable_data)  # List of normalized dfs

    #--------------------------------------------------------------------
    #  Filter data
    filter_obj = PairedSampleFilter(
        log2fold=log2fold,
        expression_min=expression_min,
        expression_max=expression_max,
        min_diff=min_diff,
        max_diff=max_diff,
        time_point_names=time_point_names,
        file_name_pairs=file_name_pairs
    )
    filtered_df_list = filter_obj.main_filter_process(split_data)
    complete_de_gene_list = filter_obj.complete_de_gene_list

    #--------------------------------------------------------------------
    # Save filter results
    # output_dir = make_default_output_dir(output_dir or None, overwrite=True)
    data_printer = DataPrinter(
        output_dir=output_dir,
        experiment_name=experiment_name,
        container_obj=container_obj,
        time_point_names=time_point_names,
        filtered_df_list=filtered_df_list,
        complete_de_gene_list=complete_de_gene_list
        )()

    bar_plotter = PairedBarPlot(
        time_point_names=time_point_names,
        output_dir=output_dir,
        experiment_name=experiment_name,
        log2fold=log2fold,
        expression_max=expression_max,
        expression_min=expression_min,
        min_diff=min_diff,
        max_diff=max_diff
    )
    bar_plotter.create_bar_plot(filter_obj.de_count_by_stage)

    scatter_plotter = ScatterPlots(
        output_dir=output_dir,
        container_obj=container_obj,
        filter_obj=filter_obj,
        experiment_name=experiment_name,
        log2fold=log2fold,
        expression_min=expression_min,
        expression_max=expression_max,
        min_diff=min_diff,
        max_diff=max_diff,
        time_point_names=time_point_names,
        scatter_min=scatter_min,
        scatter_max=scatter_max
        )
    scatter_plotter.create_scatter_plots()
    tally_plotter = TallyDe(
        output_dir=output_dir,
        container_obj=container_obj,
        experiment_name=experiment_name,
        log2fold=log2fold,
        expression_max=expression_max,
        expression_min=expression_min,
        min_diff=min_diff,
        max_diff=max_diff,
        file_name_pairs=file_name_pairs,
        time_point_names=time_point_names
    )
    tally_plotter.create_tally_plot(split_data)

    plottable_data = container_obj.normalized_df
    col_names = plottable_data.columns.tolist()

    pca_decomp = PCADecomposition(
        output_dir=output_dir,
        container_obj=container_obj,
        experiment_name=experiment_name,
        plottable_data=plottable_data,
        col_names=col_names
    )
    pca_decomp.create_pca_plot()

    #--------------------------------------------------------------------
    # Tidy up

    print("\nScript completed no errors")
    return True, ''

    # except Exception as e:
    #     print("SCRIPT FKUP")
    #     return False, e