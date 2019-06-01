

def run_splot(config_obj, correct_by_rotation):

    config_obj = config_parser(args.config)

    # load the data container_obj
    container_obj = DataContainer(config_obj)
    data, ercc_data = container_obj.parse_input()

    if args.impute:
        print('Imputation not yet implemented')

    # TODO allow this option?
    # if not args.unnorm:
    data = container_obj.normalize_file_pairs(data) # Single df of normalized data

    split_data = container_obj.split(data)  # List of normalized dfs

    # if args.svd:
    #     split_data = container_obj.remove_variance(split_data)

    if correct_by_rotation:
        split_data = container_obj.correct_via_rotation(split_data)

#--------------------------------------------------------------------
#  Filter data

    filter_obj = PairedSampleFilter(config_obj)
    filter_result = filter_obj.main_filter_process(split_data)

#--------------------------------------------------------------------
# Save filter results

    output_path = config_obj.get('data_directory', 'output')
    output_path = make_default_output_dir(output_path or None, args.overwrite)

    data_printer = DataPrinter(config_obj, container_obj or None, filter_obj or None)()

#--------------------------------------------------------------------
# Generate Plots

    print("\nPlotting data...\n")
    line_plotter = PairedDataLinePlotter(config_obj, filter_obj, data)
    fig_list = MakeFigureList(config_obj)
    line_plotter.plot_figure(figure_list=fig_list.plot_groups, plottable_data=data)

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