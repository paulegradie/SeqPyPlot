
import codecs
import os

class MakeFigureList(object):


    def __init__(self, config_obj):

        self.config = config_obj

        self.gene_list = self.input_file_parser()
        self.plot_groups = self.make_plot_groups(self.gene_list)
        assert len(self.plot_groups) > 0, "No input genes. Comon now! gene_list_parser 1"

    def input_file_parser(self):

        # Handle file encodings when you open the input file
        file_parsed = False

        input_file = self.config.get('plots', 'genelist')
        assert os.stat(input_file).st_size > 0, "Gene list for plotting was empty, gene_list_parser 2"
        
        for e in ["utf-8", "ascii", "ansi"]:
            try:
                # create list of gene names ['str1', 'str2', etc]
                input_file = codecs.open(input_file, 'r', encoding=e)
                gene_list = [row.strip().title() for row in input_file if row.strip() != '']
                file_parsed = True

            except UnicodeDecodeError:
                print("File is encoded in a format other than {}.".format(e))

            else:
                print("Parsing file using {} encoding.".format(e))
                break

        if not file_parsed:
            print("File list not parsed properly. Save file as either utf-8 or ascii encoded text.")
            raise RuntimeError
        
        return gene_list

    def make_plot_groups(self, gene_list, num_genes_per_figure=6):        
                
        # subdivide input files with more than 6 gene names.
        if len(gene_list) >= num_genes_per_figure:

            sub_list = []

            for i in range((len(gene_list) / num_genes_per_figure) + 1):
                sub_list.append(gene_list[:6])
                del (gene_list[:num_genes_per_figure])
                plot_groups = [x for x in sub_list if x != []]
        else:
            plot_groups = [[x for x in gene_list]]

        return plot_groups
