from utils.utils import make_default_dir
from csv import writer
class DataPrinter(object):

    def __init__(self, config_obj, container, analyzer):
        
        assert config_obj.get('data_directory', 'output')
        assert config_obj.get('file_names', 'prefix')

        self.output_dir = os.path.join(config_obj.get('data_directory', 'output'),
                                       config_obj.get('file_names', 'prefix'))
        self.kwargs = {'sep': '\t'}

        #self.kwargs = {'columns': cols, 'header': header?, index: sep='\t'}

    def write_ercc_data(self, ercc_df):
        
        ercc_file = "_".join([self.output_dir, "all_ercc.txt"])
        ercc_df.to_csv(ercc_file, **self.kwargs)
            
    def write_plot_data(self, plot_data_df):
        plot_data_file = "_".join([self.output_dir, 'plotter_data.txt'])
        plot_data_fd.to_csv(plot_data_file, **self.kwargs)

    def write_filtered_data(self, filtered_df_list):
        
        for stage, data in zip(self.config_obj.get('names', 'times'), filtered_df_list):
            filtered_data_file = "_".join([self.output_dir, stage, self.prefix, 'filtered_data.txt'])
            data.to_csv(filtered_data_file, **self.kwargs)

    def write_de_results(self, de_genelist):
        
        de_list_file = "_".join([self.output_dir, 'de_gene_list.txt'])
        with open(de_list_file, 'wb+') as outfile:
            result_writer = writer(outfile, delimiter='\n')
            
            for de_gene in sorted(de_genelist):
                result_writer.writerow([de_gene])