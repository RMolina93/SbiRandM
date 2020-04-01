from tkinter import *
from tkinter import ttk
from sbiRandM.sbiRandM.main import *
from sbiRandM.sbiRandM.modeller_comparison import *

window = Tk()


### MODELLER TAB ZONE

lblModPDBPath = Label(window, text="Introduce the path to the PDB file")

lblModPDBPath.grid(column=0, row=0)

lblModFastaPath = Label(window, text="Introduce the path to the Fasta file")

lblModFastaPath.grid(column=0, row=2)

lblModOutputPath = Label(window, text="Introduce the path to the Output directory")

lblModOutputPath.grid(column=0, row=4)

txtModPDBPath = Entry(window,width=20)

txtModPDBPath.grid(column=1, row=0)

txtModFastaPath = Entry(window,width=20)

txtModFastaPath.grid(column=1, row=2)

txtModOutputPath = Entry(window,width=20)

txtModOutputPath.grid(column=1, row=4)

var = IntVar()

radioMod = Radiobutton(window, text="Use Modeller", variable=var, value=1)
radioMod.grid(column=0, row=6)

radioSup = Radiobutton(window, text='Use Superimp', variable=var, value=2)
radioSup.grid(column=1, row=6)

### FUNCTION ZONE
def clicked():
   #Codigo de ejecucion modeller
   args = dict()
   args['folder'] = txtModPDBPath.get()
   args['fasta_seq'] = txtModFastaPath.get()
   args['output_folder'] = txtModOutputPath.get()
   
   print(args)
   if var.get() == 2:
      #Por superimposicion
      mainSuperimp(args)
   elif var.get() == 1:
      mainMod(args)
### END OF FUNCTION ZONE

btnMod = Button(window, text='Run!', command=clicked)

btnMod.grid(column=3, row=6)

window.mainloop()






