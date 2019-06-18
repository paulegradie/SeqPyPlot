from main_app.seqpyplot.parsers.gene_list_parser import MakeFigureList
from main_app.seqpyplot.plot.paired_line_plotter import PairedDataLinePlotter
from seqpyplot.parsers.config_parser import config_parser
from seqpyplot.container.data_container import DataContainer
from seqpyplot.analyzer.paired_sample_filter import PairedSampleFilter

from pathlib import Path
import os

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

# def write_textbox_genelist_to_tmp(genelist, tmp_plot_dir):
#     output_path = os.path.join(tmp_plot_dir, 'genelist.txt')
#     with open(output_path, 'w') as fout:
#         for gene in genelist:
#             fout.write(gene + '\n')
#     return output_path


def generate_plots(
    request,
    data_type,
    tmp_output_dir,
    min_diff,
    max_diff
):

    num_conditions = check_num_series(request)
    processed_gene_list = collect_gene_list(request)
    plottable_data_path = upload_plottable_data(request, tmp_output_dir)
    complete_de_genelist = upload_de_genelist(request, tmp_output_dir)

    # the params need t come from the config
    config_path = collect_plot_config_path(request, tmp_output_dir)
    config_obj = config_parser(config_path)
    time_point_names = config_obj.get('names', 'times')
    control_sample_names = config_obj.get('names', 'controls')
    treated_sample_names = config_obj.get('names', 'treated')
    log2fold = config_obj.getfloat('params', 'log2fold')
    expression_min = config_obj.get('params', 'expression_min')
    min_diff =  config_obj.getint('params', 'min_diff')
    max_diff =  config_obj.getint('params', 'max_diff')
    experiment_name = config_obj.get('names', 'experiment_name')


    sample_names = control_sample_names + treated_sample_names
    file_name_pairs = [(x, y) for x, y in zip(control_sample_names, treated_sample_names)]
    num_file_pairs = len(file_name_pairs)

    container_obj = DataContainer(
        data_directory=None,  # not used here
        data_paths=None,  # not used here
        sample_names=sample_names,
        time_point_names=time_point_names,
        file_name_pairs=file_name_pairs,
        num_file_pairs=num_file_pairs,
        num_components=None,  # this is for svd variance reduction. Not necessary atm - not implemented
        data_type=data_type)

    normalized_data, ercc_data = container_obj.parse_plotter_data(plottable_data_path)

    # import pdb; pdb.set_trace()
    print("\nPlotting data...\n")
    line_plotter = PairedDataLinePlotter(
        normalized_df=normalized_data,
        complete_de_gene_list=complete_de_genelist,
        log2fold=log2fold,
        expression_min=expression_min,
        time_point_names=time_point_names,
        condition_labels=['Control'] if num_conditions == 1 else ['Control', 'Treated'], # only 1 or two supported right now.
        experiment_name=experiment_name,
        min_diff = min_diff,
        max_diff = max_diff
    )

    fig_list = MakeFigureList(genelist=processed_gene_list)
    line_plotter.plot_figure(figure_list=fig_list.plot_groups, plottable_data=normalized_data)
    image_files = [os.path.join(tmp_output_dir, x) for x in os.listdir(tmp_output_dir) if x.endswith('png')]
    return image_files
