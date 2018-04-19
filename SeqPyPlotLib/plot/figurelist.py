
import codecs

class MakeFigureList(object):
    """
    Function to parse the input gene list.
    """

    def __init__(self, args):
        """
        :param args: Args object
        """
        self.args = args
        self.gene_list = []  # The full gene list from the input file
        self.figure_list, self.gene_list = self.input_list_parser()  # nested list - gene_list broken in to gruops of 6

    def input_list_parser(self):

        """
        # type: (input_file object)
        :return: 2D nested list object

        """
        # Handle file encodings when you open the input file
        for e in ["utf-8", "ascii", "ansi"]:
            try:
                # create list of gene names ['str1', 'str2', etc]
                genes_in = codecs.open(self.args.gene_list, 'r', encoding=e)
                for row in genes_in:
                    if row.rstrip() != '':
                        self.gene_list.append(str(row.rstrip().capitalize()))
            except UnicodeDecodeError:
                print("File is encoded in a format other than {}.".format(e))

            else:
                print("Parsing file using {} encoding.".format(e))
                break

        # Handle the case where ultimately the file can't be opened
        if len(self.gene_list) == 0:
            print("File list not parsed properly. Save file as either utf-8 or ascii encoded text.")
            raise RuntimeError

        # subdivide in put files with more than 6 gene names.
        if len(self.gene_list) >= 6:

            sub_list = []

            for i in range((len(self.gene_list) / 6) + 1):
                sub_list.append(self.gene_list[:6])
                del (self.gene_list[:6])
                self.figure_list = [x for x in sub_list if x != []]

            return self.figure_list, self.gene_list

        else:
            self.figure_list = [[x for x in self.gene_list]]

            return self.figure_list, self.gene_list

    def figure_list_length(self):
        if len(self.figure_list) == 0:
            print("Your gene list is currently empty or couldn't be opened.")
        else:
            return len(self.figure_list)