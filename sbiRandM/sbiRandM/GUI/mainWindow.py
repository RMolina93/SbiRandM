from tkinter import *
from tkinter import ttk
from sbiRandM.sbiRandM.superimp import *
from sbiRandM.sbiRandM.modeller_comparison import *
from tkinter import filedialog

class ventana():

   
    def __init__(self , master):


        ### MODELLER TAB ZONE
      
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

   ### FUNCTION ZONE
    def clicked(self):
        #Codigo de ejecucion modeller
        args = dict()
        args['folder'] = self.txtModPDBPath.get()
        args['fasta_seq'] = self.txtModFastaPath.get()
        args['output_folder'] = self.txtModOutputPath.get()
        args['verbose'] = self.verb.get()
      
        if self.var.get() == 2:
           #Por superimposicion
           mainSuperimp(args)
        elif self.var.get() == 1:
           mainMod(args)

        

    ### END OF FUNCTION ZONE

    def browsePDB(self):
        
        directory = filedialog.askdirectory()
        self.txtModPDBPath.delete(0 , 'end')
        self.txtModPDBPath.insert(0,directory)

    def browseFasta(self):
   
       file = filedialog.askopenfilename()
       self.txtModFastaPath.delete(0 , 'end')
       self.txtModFastaPath.insert(0,file )

    def browseOutput(self):
   
       directory = filedialog.askdirectory()
       self.txtModOutputPath.delete(0 , 'end')
       self.txtModOutputPath.insert(0,directory)
       
if __name__ == "__main__":
    root = Tk()
    
    ss = ventana(root)
    
    root.title("sbiRandM complex builder")
    root['bg'] = "#EAEAEA"
    root.resizable(False, False)
    root.mainloop()







