import numpy as np
from tkinter import *
from tkinter import ttk
import matplotlib as plt

R = 6350000	# Radius of NVA-3

def initial_altitude_changed(*args):
	print("initial altitude was changed")
	delta_longitude.set(180 * (1 - np.power((2*R + (a_i.get() + a_f.get())*1000)/(2*R+2*a_f.get()*1000), 3/2)))
	print(f"delta longitude is now {delta_longitude.get()}")
	return
	
def final_altitude_changed(*args):
	print("final altitude was changed")
	delta_longitude.set(180 * (1 - np.power((2*R + (a_i.get() + a_f.get())*1000)/(2*R+2*a_f.get()*1000), 3/2)))
	print(f"delta longitude is now {delta_longitude.get()}")
	return

def text_changed(*args):
	print("E")
	return

root = Tk()

a_i = DoubleVar(value = 185.0)
a_f = DoubleVar(value = 486.0)
delta_longitude = DoubleVar()
delta_longitude.set(180 * (1 - np.power((2*R + (a_i.get() + a_f.get())*1000)/(2*R+2*a_f.get()*1000), 3/2)))

frm = ttk.Frame(root, padding=10)
frm.grid()
ttk.Label(frm, text="Initial altitude").grid(column=0, row=0)
a_i_Entry = ttk.Entry(frm, textvariable = a_i)
a_i_Entry.grid(column = 1, row = 0)
a_i.trace_add("write", initial_altitude_changed)
ttk.Label(frm, text="km").grid(column = 2, row = 0)

ttk.Label(frm, text="Final altitude").grid(column=0, row=1)
a_f_Entry = ttk.Entry(frm, textvariable = a_f)
a_f_Entry.grid(column = 1, row = 1)
a_f.trace_add("write", final_altitude_changed)
ttk.Label(frm, text="km").grid(column = 2, row = 1)

ttk.Label(frm, text = "Delta Longitude").grid(column = 0, row = 2)
delta_longitude_label = ttk.Label(frm, textvariable = delta_longitude)
delta_longitude_label.grid(column = 1, row = 2)
ttk.Label(frm, text="⁰").grid(column = 2, row = 2)

root.mainloop()

# ~ a_i = float(input("Enter initial altitude in km: ")) *1000
# ~ a_f = float(input("Enter final altitude in km: ")) * 1000

# ~ delta_longitude = 180 * (1 - np.power((2*R + a_i + a_f)/(2*R+2*a_f), 3/2))

# ~ print(f"Burn at relative longitude of {delta_longitude}")
