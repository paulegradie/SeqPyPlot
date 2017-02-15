from tkinter import *
from tkFileDialog import askopenfilename
import Tkconstants

class Interface:
    def __init__(self, master):

        self.master = master
        master.title("SeqPyPlot v0.2 GUI")
        master.geometry('680x500')

        self.button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}

        # main script string
        self.callstring = dict()

        # defining options for opening a directory
        self.dir_opt = options = {}
        options['initialdir'] = 'C:\\'
        options['mustexist'] = False
        options['parent'] = root
        options['title'] = 'This is a title'

    def set_entry(self, labels):

        row = 0
        col = 2
        for label, default_name, call, default in labels:
            optionlabel = Label(self.master, text=label)
            optionlabel.grid(row=row, column=col)
            row += 1
            e = Entry(self.master)
            e.insert(END, default_name)
            e.grid(row=row, column=col)

            row += 1

    def file_handler(self, opt, vars, pos, row, col):
        filepath = askopenfilename()
        self.callstring[opt] = filepath
        Label(self.master, text='...' + filepath[-20:]).grid(row=row, column=col+1)
        var[pos].get()

    def set_file_loader(self, loader_list):

        row = 0
        col = 3
        buttons = []
        vars = []
        pos = 0

        for name, opt, default in loader_list:
            Label(self.master, text=name).grid(row=row, column=col)
            row += 1

            vars.append(StringVar(self.master))
            vars[pos].set(default)

            buttons.append(Button(self.master, text=name, command=lambda: self.file_handler(opt, vars, pos, row, col), padx=10, pady=5))
            buttons[pos].grid(row=row, column=col, padx=10)

            row += 1
            pos += 1

    def set_dropdown(self, menu_list):
        col = 1
        row = 0
        for label, options, default_name, call, default in menu_list:

            menulabel = Label(self.master, text=label)
            menulabel.grid(row=row, column=col)
            row += 1

            variable = StringVar(self.master)
            variable.set(default)
            menu = OptionMenu(self.master, variable, *options)
            menu.grid(row=row, column=col)
            row += 1


if __name__ == "__main__":

    root = Tk()
    seqPyPlot_GUI = Interface(root)
    entry_list = [('Stages', 'None', '-time ', 'None'),
                  ('Output folder name', 'None', '-out ', 'None'),
                  ('Conditions', 'Series1,Series2', '-condition ', 'Series1,Series2'),
                  ('Lower Expr Limit', '10', '-low ', '10'),
                  ('Upper Expr Limit', '2000', '-hi ', '2000'),
                  ('Min Expr Diff', '3.5', '-dif ', '3.5'),
                  ('Log2fold Cutoff', '0.7', '-log2 ', '0.7'),
                  ('File Prefix', 'SeqPyPlot_out', '-prefix ', 'SeqPyPlot_out')]

    # Label, options, default_name, call, default
    menu_list = (('Data Type', ('htseq', 'cuffnorm'), 'htseq', '-data_type ', 'htseq'),
                 ('Number of Plots', ('1', '2'), '2', '-num ', '2'),
                 ('Remove Transient/Off Genes', ('Yes', 'No'), 'No', '-remove ', 'False'),
                 ('Tally flagged Genes', ('Yes', 'No'), 'Yes', '-tally ', 'True'),
                 ('Write out Filter and Plot Files', ('Yes', 'No'), 'Yes', '-report ', 'True'),
                 ('Write out ERCC data', ('Yes', 'No'), 'No', '-ercc ', 'False'))

    loader_list = [('Filter Results', '-fr ', 'None'),
                   ('Raw Data', '-raw_data ', 'None'),
                   ('Plot Data', '-plot_data ', 'None'),
                   ('Gene List to Plot', '-gene_list ', 'None'),
                   ('Custom DE List', '-results ', 'None')]

    seqPyPlot_GUI.set_entry(entry_list)
    seqPyPlot_GUI.set_dropdown(menu_list)
    seqPyPlot_GUI.set_file_loader(loader_list)
    Button(root, text="TEST").grid(row=10, column=10)

    print seqPyPlot_GUI.callstring

    root.mainloop()



