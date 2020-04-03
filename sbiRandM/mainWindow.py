from tkinter import *
from tkinter import ttk
import os
from tkinter import messagebox
from .superimp import *

try:
    from .modeller_comparison import *
except:
    pass

from .exceptions import *
from tkinter import filedialog
import logging

os.environ['TK_SILENCE_DEPRECATION'] = "1"

class window():
    """
    This class is defined as the main window of the Graphical User Interface.
    """
   
    def __init__(self , master):

        logging.basicConfig(filename='example.log', level=logging.DEBUG)
      
        self.lblModPDBPath = ttk.Label(master, text="Introduce the path to the PDB directory")
        self.lblModPDBPath.grid(column=0, row=0)
      
        self.lblModFastaPath = ttk.Label(master, text="Introduce the path to the fasta file")
        self.lblModFastaPath.grid(column=0, row=2)
      
        self.lblModOutputPath = ttk.Label(master, text="Introduce the path to the output directory")
        self.lblModOutputPath.grid(column=0, row=4)
      
        self.txtModPDBPath = ttk.Entry(master,width=40)
        self.txtModPDBPath.grid(column=1, row=0)
      
        self.txtModFastaPath = ttk.Entry(master,width=40)
        self.txtModFastaPath.grid(column=1, row=2)
      
        self.txtModOutputPath = ttk.Entry(master,width=40)
        self.txtModOutputPath.grid(column=1, row=4)
        
        self.lblOut = ttk.Label(master , text="You can follow the results in the Terminal!")
        self.lblOut.grid(column=1 , row=9)
        
        self.var = IntVar()
        self.radioMod = ttk.Radiobutton(master , text="Use Modeller" , variable=self.var , value=1)
        self.radioMod.grid(column=1 , row=6)
        self.radioSup = ttk.Radiobutton(master , text='Use Superimposition' , variable=self.var , value=2)
        self.radioSup.grid(column=0 , row=6)
        self.var.set(2)

        self.verbose = BooleanVar()
        self.verb = ttk.Checkbutton(master, text="verbose", variable=self.verbose)
        self.verb.grid(column=2, row=6)

        self.btnMod = ttk.Button(master , text='Run!' , command=self.clicked, width=20)
        self.btnMod.grid(column=3 , row=6)

        self.btn1 = ttk.Button(master , text='Browse PDB' , command=self.browsePDB , width=10)
        self.btn1.grid(column=3 , row=0)

        self.btn2 = ttk.Button(master , text='Browse fasta' , command=self.browseFasta , width=10)
        self.btn2.grid(column=3 , row=2)

        self.btn3 = ttk.Button(master , text='Select output' , command=self.browseOutput , width=10)
        self.btn3.grid(column=3 , row=4)


    def clicked(self):
        args = dict()
        if (self.txtModPDBPath.get() == "" or self.txtModFastaPath.get() == "" or self.txtModOutputPath.get() == ""):
           messagebox.showinfo("Error (But don't worry!)" ,
                               "You didn't enter a path or file for the arguments! Please feel free to try again.")
           self.clear()

        else:
   
           args['folder'] = self.txtModPDBPath.get()
           args['fasta_seq'] = self.txtModFastaPath.get()
           args['output_folder'] = self.txtModOutputPath.get()
           args['verbose'] = self.verbose.get()
           try:
                if self.var.get() == 2:
                   mainSuperimp(args)   #Call superimposition-based approach
                elif self.var.get() == 1:
                   mainMod(args)        #Call Modeling-based approach

           except FileNotFoundError:
                messagebox.showinfo("Error (But don't worry!)" , "It seems like one of the paths that you introduced is not valid. Please try again!")
                logging.error('It seems like one of the paths that you introduced is not valid. Please try again!')
                self.clear()
           except IsADirectoryError:
                messagebox.showinfo("Error (But don't worry!)" , "It seems like one of the paths that you introduced is a directory instead of a file. Please try again!")
                logging.error('It seems like one of the paths that you introduced is a directory instead of a file. Please try again!')
                self.clear()
           except FilesDontMatchException:
                messagebox.showinfo("Error (But don't worry!)", "The files you introduced are not correct")
                logging.error('The files you introduced are not correct')
                self.clear()
           except FastaRaroException:
                messagebox.showinfo("Error (But don't worry!)" , "The fasta file does not match with the pdbs of the protein!")
                logging.error('The fasta file does not match with the pdbs of the protein!')
                self.clear()
           except BadFastaException:
                messagebox.showinfo("Error (But don't worry!)" ,
                                    "The fasta file is not formatted or configured correctly")
                logging.error('The fasta file does not match with the pdbs of the protein!')
                self.clear()

    def clear(self):
        """
        This function clears the paths in the GUI interface
        """
        self.txtModPDBPath.delete(0 , 'end')
        self.txtModFastaPath.delete(0 , 'end')
        self.txtModOutputPath.delete(0 , 'end')


    def browsePDB(self):
        """
        This function makes a System Call to let the user
        browse for a folder. 
        Then it inserts the directory into the GUI field.
        """
        
        directory = filedialog.askdirectory()
        self.txtModPDBPath.delete(0 , 'end')
        self.txtModPDBPath.insert(0,directory)

    def browseFasta(self):
        """
        This function makes a System Call to let the user
        browse for a fasta file. 
        Then it inserts the file into the GUI field.
        """

        file = filedialog.askopenfilename()
        self.txtModFastaPath.delete(0 , 'end')
        self.txtModFastaPath.insert(0,file )

    def browseOutput(self):
        """
        This function makes a System Call to let the user
        browse for an output directory. 
        Then it inserts the file into the GUI field.
        """
        directory = filedialog.askdirectory()
        self.txtModOutputPath.delete(0 , 'end')
        self.txtModOutputPath.insert(0,directory)

    def validate_input(self, input):
        """
        This function validates that the inputs
        in the GUI fields are strings.
        """
        if type(input) == str:
            return True
        return False


if __name__ == "__main__":

    print ("Welcome to the GUI Interface of SbiRandM Complex Builder. Please fill the options to proceed.")
    root = Tk()
    ss = window(root)
    root.title("sbiRandM Complex Builder")
    root['bg'] = "#EAEAEA"
    root.resizable(False, False)
    root.mainloop()

def exe():
   print("Welcome to the GUI Interface of SbiRandM Complex Builder. Please fill the options to proceed.")
   root = Tk()
   ss = window(root)
   root.title("sbiRandM Complex Builder")
   root['bg'] = "#EAEAEA"
   root.resizable(False , False)
   root.mainloop()




