
from Tkinter import *
from tkFileDialog import askopenfile
from tkFileDialog import asksaveasfilename
from tkFileDialog import *

#define GUI funcitons for opening presets, and saving presets
def new_preset():
    print "New File!"


def open_preset():
    name = askopenfilename()

def save_preset():
    pass

def About():
    about_message = "SeqPyPlot usage and Docs can be found at\nhttps://github.com/paulgradie/SeqPyPlot"
    msg = Message(root, text=Welcome)
    msg.config(bg='white', font=('arial', 16, 'italic'))
    msg.pack(width=20)


root = Tk()
root.geometry("720x500+50+50")

menu = Menu(root)
root.config(menu=menu)
filemenu = Menu(menu)
menu.add_cascade(label="Presets", menu=filemenu)
filemenu.add_command(label="Create", command=new_preset())
filemenu.add_command(label="Open...", command=open_preset())
filemenu.add_separator()
filemenu.add_command(label="Exit", command=root.quit)

helpmenu = Menu(menu)
menu.add_cascade(label="Help", menu=helpmenu)
helpmenu.add_command(label="About...", command=About)
x = StringVar()

mainloop()


sys.exit()
from tkFileDialog import askopenfilename





menu = Menu(root)
root.config(menu=menu)
filemenu = Menu(menu)
menu.add_cascade(label="File", menu=filemenu)
filemenu.add_command(label="New", command=NewFile)
filemenu.add_command(label="Open...", command=OpenFile)
filemenu.add_separator()
filemenu.add_command(label="Exit", command=root.quit)

helpmenu = Menu(menu)
menu.add_cascade(label="Help", menu=helpmenu)
helpmenu.add_command(label="About...", command=About)


Label(root, text="SeqPyPlot", font=('arial', 25, 'bold'), bg="white", fg='black', ).pack()


Welcome= "Welcome to SeqPyPlt simple analysis tool for Time Series RNAseq."
msg = Message(root, text = Welcome)
msg.config(bg='white', font=('arial', 16, 'italic'))
msg.pack()

def show_entry_fields():
   print("First Name: %s\nLast Name: %s" % (e1.get(), e2.get()))

Label(root, text="First Name").grid(row=0)
Label(root, text="Last Name").grid(row=1)

e1 = Entry(root)
e2 = Entry(root)

e1.grid(row=0, column=1)
e2.grid(row=1, column=1)

Button(root, text='Quit', command=root.quit).grid(row=3, column=0, sticky=W, pady=4)
Button(root, text='Show', command=show_entry_fields).grid(row=3, column=1, sticky=W, pady=4)

def show_entry_fields():
   print("First Name: %s\nLast Name: %s" % (e1.get(), e2.get()))

Label(root, text="First Name").grid(row=0)
Label(root, text="Last Name").grid(row=1)

e1 = Entry(root)
e2 = Entry(root)

e1.grid(row=0, column=1)
e2.grid(row=1, column=1)

Button(root, text='Quit', command=root.quit).grid(row=3, column=0, sticky=W, pady=4)
Button(root, text='Show', command=show_entry_fields).grid(row=3, column=1, sticky=W, pady=4)

mainloop()