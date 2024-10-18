import numpy as np
from tkinter import *
from tkinter import ttk

import matplotlib.pyplot as plt
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure
from matplotlib.patches import Circle


R = 6350000	# Radius of NVA-3
mu = 398600441800000 # Standard gravitational parameter of NVA-3 (assuming identical mass to earth)

def recalc_burn_params(*args):
	try:
		r_i = a_i.get()*1000 + R
		r_f = a_f.get()*1000 + R
		delta_longitude_degrees.set(180 * (1 - np.power((r_i + r_f)/(2*r_f), 3/2)))	# the formula :O
		delta_longitude_degrees.set(round(delta_longitude_degrees.get() % (360*np.sign(delta_longitude_degrees.get())), 2))	# can be optimised
		if abs(delta_longitude_degrees.get()) > 180:
			delta_longitude_degrees.set(360 - abs(delta_longitude_degrees.get()))	# can DEFINITELY be optimised
		
		burn_distance.set(round(np.sqrt(r_i**2+r_f**2 - 2*r_i*r_f*np.cos(np.radians(delta_longitude_degrees.get())))/1000, 2))
		
		transfer_duration.set(round((np.pi*np.sqrt((((r_i+r_f)/2)**3)/mu))/60, 1))
		
		c = -(r_f-r_i)/2
		
		a = (r_i+r_f)/2
		
		update_plot(r_f, r_i, a, c)
		
	except TclError:
		return

def update_plot(r_f, r_i, a, c):
	theta = np.linspace(0, 2*np.pi, num = 1000)
	r_transfer = (-a**2 +c**2)/(-a + c*np.cos(theta))
	
	initial_orbit.set_data(theta, np.full_like(theta, r_i))
	final_orbit.set_data(theta, np.full_like(theta, r_f))
	transfer_orbit.set_data(theta[:theta.size//2], r_transfer[:theta.size//2])
	transfer_orbit_trace.set_data(theta[theta.size//2:], r_transfer[theta.size//2:])
	ship.set_data(0, r_i)
	station.set_data(np.radians(delta_longitude_degrees.get()), r_f)
	rendezvous_point.set_data(np.pi, r_f)
	
	ax.set_rmax(max(r_i, r_f)*1.1)
	
	canvas.draw()
	return

def altitude_changed(*args):
	recalc_burn_params()
	return

root = Tk()
root.wm_title("Inter-orbital Rendezvous Calculator")

a_i = DoubleVar(value = 185.0)
a_f = DoubleVar(value = 1058.0)
delta_longitude_degrees = DoubleVar()
burn_distance = DoubleVar(value = 1000)
transfer_duration = DoubleVar(value = 1)

frm = ttk.Frame(root, padding=10)
frm.grid()

ttk.Label(frm, text="Initial altitude").grid(column=0, row=0)
a_i_Entry = ttk.Entry(frm, textvariable = a_i)
a_i_Entry.grid(column = 1, row = 0)
a_i.trace_add("write", altitude_changed)
ttk.Label(frm, text="km").grid(column = 2, row = 0)

ttk.Label(frm, text="Final altitude").grid(column=0, row=1)
a_f_Entry = ttk.Entry(frm, textvariable = a_f)
a_f_Entry.grid(column = 1, row = 1)
a_f.trace_add("write", altitude_changed)
ttk.Label(frm, text="km").grid(column = 2, row = 1)

ttk.Label(frm, text = "Delta Longitude").grid(column = 0, row = 2)
delta_longitude_degrees_label = ttk.Label(frm, textvariable = delta_longitude_degrees)
delta_longitude_degrees_label.grid(column = 1, row = 2)
ttk.Label(frm, text="⁰").grid(column = 2, row = 2)

ttk.Label(frm, text = "Distance to target").grid(column = 0, row = 3)
burn_distance_label = ttk.Label(frm, textvariable = burn_distance)
burn_distance_label.grid(column = 1, row = 3)
ttk.Label(frm, text="km").grid(column = 2, row = 3)

ttk.Label(frm, text = "Transfer duration").grid(column = 0, row = 4)
transfer_duration_label = ttk.Label(frm, textvariable = transfer_duration)
transfer_duration_label.grid(column = 1, row = 4)
ttk.Label(frm, text="min").grid(column = 2, row = 4)

# ~ set up canvas for orbit plots

theta = np.linspace(0, 2*np.pi, num = 1000)
r_i = a_i.get()*1000 + R
r_f = a_f.get()*1000 + R

# ~ elliptical orbit parameters

c = -(r_f-r_i)/2

a = (r_i+r_f)/2

r_transfer = (-a**2 +c**2)/(-a + c*np.cos(theta))

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize = [9, 9], layout = "tight")
fig.set_facecolor("midnightblue")
fig.set_frameon(True)
fig.set_edgecolor("darkgrey")
ax.spines['polar'].set_visible(False)
ax.set_facecolor("midnightblue")
planet = Circle((0,0), R, transform=ax.transData._b, facecolor = "teal", edgecolor = "dimgrey")
ax.add_artist(planet)
initial_orbit, = ax.plot(theta, np.full_like(theta, r_i), color = "red")
final_orbit, = ax.plot(theta, np.full_like(theta, r_f), color = "skyblue")
transfer_orbit, = ax.plot(theta[:theta.size//2], r_transfer[:theta.size//2], color = "limegreen")
transfer_orbit_trace, = ax.plot(theta[theta.size//2:], r_transfer[theta.size//2:], color = "limegreen", linestyle = 'dashed')
ship, = ax.plot(0, r_i, color = "white", marker = "^", markersize = 15)
station, = ax.plot(np.radians(delta_longitude_degrees.get()), r_f, color = "skyblue", marker = "D", markersize = 12)
rendezvous_point, = ax.plot(np.pi, r_f, color = "orange", marker = "X", markersize = 12)

ax.set_rmax(max(r_i, r_f)*1.1)
ax.set_rticks([])	# No radial ticks
ax.set_rgrids([])	# No radial grids
ax.set_theta_offset(np.pi)
ax.set_theta_direction(-1)
ax.set_thetagrids([])
ax.grid(False)

canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.

recalc_burn_params()

canvas.get_tk_widget().grid(column = 4, row = 0, padx = 4, pady = 4)


def on_closing():
	root.destroy()
	quit()

root.protocol("WM_DELETE_WINDOW", on_closing)
root.mainloop()

# ~ a_i = float(input("Enter initial altitude in km: ")) *1000
# ~ a_f = float(input("Enter final altitude in km: ")) * 1000

# ~ delta_longitude_degrees = 180 * (1 - np.power((2*R + a_i + a_f)/(2*R+2*a_f), 3/2))

# ~ print(f"Burn at relative longitude of {delta_longitude_degrees}")
