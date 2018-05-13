
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
    
    # load the data container
    dc = DataContainer(config)
    raw_df, ercc_data = dc._parse_input_()
    
    normalized_df = dc.reorder_cols(dc.normalize_file_pairs())
    complete_gene_list = normalized_df.index.tolist()

    split_normalized_dfs = dc.split()

    if dc.args.impute_by_nieghbors:
        normalized_df = dc._average_flanking_()

    if write_csv:
        write_to_csv(raw_df, 'test_raw_df.txt')
        write_to_csv(raw_df, 'test_ercc_df.txt')
        for df in split_normalized_dfs:
            write_to_csv(df, 'test_raw_normed_df.txt')