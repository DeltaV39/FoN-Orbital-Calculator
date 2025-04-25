from tkinter import *
from tkinter import ttk
from hohmann import HohmannCalculator
from discrete import DiscreteCalculator
from reentry import ReentryCalculator

root = Tk()
root.wm_title("Inter-orbital Rendezvous Calculator")
notebook = ttk.Notebook(root)

tab_hohmann = HohmannCalculator(notebook)
tab_discrete = DiscreteCalculator(notebook)
tab_reentry = ReentryCalculator(notebook)

notebook.add(tab_hohmann, text="Hohmann Transfer")
notebook.add(tab_discrete, text="Custom Transfer")
notebook.add(tab_reentry, text="Reentry")

notebook.pack(expand=True, fill='both')

def on_closing():
    root.destroy()
    quit()

root.protocol("WM_DELETE_WINDOW", on_closing)
root.mainloop()
