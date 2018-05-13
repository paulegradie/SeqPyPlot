
from SeqPyPlot.seqpyplot.container.data_container import DataContainer
from argparse import ArgumentParser

if __name__ == "__main__":
    """
    Script to normalize data frame independantly
    """
    usage = """
    This is how you use this program
    """
    name = "SeqpyPlot v0.4"
    parser = ArgumentParser(usage=usage, description=name)
    
    parser.add_argument('-c', '--config', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, default='default')
    parser.add_argument('-i', '--impute', type=bool, action='store_true')
    args = parser.parse_args()


    # load the data container
    dc = DataContainer(args.config)
    raw_df, ercc_data = dc._parse_input_()
    
    normalized_df = dc.reorder_cols(dc.normalize_file_pairs())
    complete_gene_list = normalized_df.index.tolist()

    split_normalized_dfs = dc.split()

    if dc.args.impute_by_nieghbors:
        normalized_df = dc._average_flanking_()

    if write_csv:
        write_to_csv(raw_df, '{}_raw_df.txt'.format(args.output))
        write_to_csv(ercc_data, '{}_ercc_df.txt'.format(args.output))
        for df in split_normalized_dfs:
            write_to_csv(df, '{}_raw_normed_df.txt'.format(args.output))