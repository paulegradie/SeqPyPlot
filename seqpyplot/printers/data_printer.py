import os
from csv import writer

from seqpyplot.utils.utils import make_default_output_dir

import os
from csv import writer
from datetime import datetime

TIME_STAMP = str(datetime.now()).split()[0]


# def write_sing_col_list(data_list, output_path, filename):
#     """
#     Write a list to a file as a single column
#     Intended for a gene list
#     """
#     output_path = os.path.join(output_path, "_".join([filename, TIMESTAMP]))
#     with open(output_path, 'wb+') as outfile:
#         result_writer = writer(outfile, delimiter='\n')
        
#         for de_gene in sorted(data_list):
#             result_writer.writerow([de_gene])




# def write_df_to_csv(data, file_path, suffix=None):
#     """ write dataframe to csv """

#     assert str(type(data)) == "<class 'pandas.core.frame.DataFrame'>"
#     file_name = os.path.join(file_path, suffix if suffix else '')
#     data.to_csv(file_name, sep='\t')
#     return None




class DataPrinter(object):

    def __init__(self, config_obj, container, filter_obj):

        assert config_obj.get('data_directory', 'output')
        assert config_obj.get('file_names', 'prefix')

        self.output_dir = config_obj.get('data_directory', 'output')
        self.prefix = config_obj.get('names', 'experiment_name')

        self.config_obj = config_obj
        self.container = container
        self.filter_obj = filter_obj

        self.kwargs = {'sep': '\t'}

    def __call__(self):

        foos_to_print = [
            self.write_ercc_data,
            self.write_normalized_data,
            self.write_filtered_data,
            self.write_complete_de_list
            ]
        
        for foo in foos_to_print:
            try:
                foo()
            except NameError:
                print("No {} data to write".format(foo.__name__))
                pass

    def write_ercc_data(self):

        ercc_data = self.container.ercc_df

        output_name = '_'.join([self.prefix, "all_ercc.spp"])
        output_path = os.path.join(self.output_dir, output_name)

        ercc_data.to_csv(output_path, **self.kwargs)
            
    def write_normalized_data(self):

        normalized_data = self.container.normalized_data

        output_name = '_'.join([self.prefix, "normalized_data.spp"])
        output_path = os.path.join(self.output_dir, output_name)

        normalized_data.to_csv(output_path, **self.kwargs)

    def write_filtered_data(self):
        
        filtered_df_list = self.filter_obj.filtered_df_list

        for stage, data in zip(self.config_obj.getlist('names', 'times'), filtered_df_list):
            
            output_name = '_'.join([self.prefix, stage, 'filtered_data.spp'])
            output_path = os.path.join(self.output_dir, output_name)            
            data.to_csv(output_path, **self.kwargs)

    def write_complete_de_list(self):

        de_gene_list = list(self.filter_obj.complete_de_gene_list)
        output_name = "_".join([self.prefix, 'de_gene_list.spp'])
        output_path = os.path.join(self.output_dir, output_name)            

        with open(output_path, 'w+') as outfile:
            result_writer = writer(outfile, delimiter=' ')
            
            for de_gene in sorted(de_gene_list):
                result_writer.writerow([de_gene])
