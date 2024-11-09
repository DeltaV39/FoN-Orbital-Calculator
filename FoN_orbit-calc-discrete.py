import numpy as np
from tkinter import *
from tkinter import ttk
import scipy.constants

import matplotlib.pyplot as plt
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

R = 6350	# Radius of NVA-3 in km
M = 5.9722e24
mu = 398600441800000 # Standard gravitational parameter of NVA-3 (assuming identical mass to earth)

def altitude_changed(*args):
	try:
		r_i = a_i.get() + R
		r_f = a_f.get() + R
	
	except TclError:
		return
	
	# update plot with new altitudes, rerun targeting algo if implemented
	return

# ~ def Hohmann_transfer(r_i, r_f):
	# ~ delta_longitude_degrees.set(180 * (1 - np.power((r_i + r_f)/(2*r_f), 3/2)))	# the formula :O
	# ~ delta_longitude_degrees.set(round(delta_longitude_degrees.get() % (360*np.sign(delta_longitude_degrees.get())), 2))	# can be optimised
	# ~ if abs(delta_longitude_degrees.get()) > 180:
		# ~ delta_longitude_degrees.set(360 - abs(delta_longitude_degrees.get()))	# can DEFINITELY be optimised
	
	# ~ burn_distance.set(round(np.sqrt(r_i**2+r_f**2 - 2*r_i*r_f*np.cos(np.radians(delta_longitude_degrees.get())))/1000, 2))
	
	# ~ transfer_duration.set(round((np.pi*np.sqrt((((r_i+r_f)/2)**3)/mu))/60, 1))
	
	# ~ c = -(r_f-r_i)/2
	
	# ~ a = (r_i+r_f)/2
	
	# ~ r_transfer = (-a**2 +c**2)/(-a + c*np.cos(theta))
	
	# ~ initial_orbit.set(radius = r_i)
	# ~ final_orbit.set(radius = r_f)
	# ~ transfer_orbit.set_data(theta[:theta.size//2], r_transfer[:theta.size//2])
	# ~ transfer_orbit_trace.set_data(theta[theta.size//2:], r_transfer[theta.size//2:])
	# ~ ship.set_data(0, r_i)
	# ~ station.set_data(np.radians(delta_longitude_degrees.get()), r_f)
	# ~ rendezvous_point.set_data(np.pi, r_f)
	# ~ ax.set_rmax(max(r_i, r_f)*1.1)
	
	# ~ canvas.draw()
	# ~ return

def ang_speed_from_alt(altitude):
	
	speed = np.sqrt((scipy.constants.G*M)/numpy.pow(altitude,3))
	
	return speed

root = Tk()
root.wm_title("Inter-orbital Rendezvous Calculator")

a_i = DoubleVar(value = 185.0)
a_f = DoubleVar(value = 1058.0)
r_i = 185.0 + R
r_f = 1058.0 + R
delta_longitude_degrees = DoubleVar()
burn_distance = DoubleVar(value = 1000)
transfer_duration = DoubleVar(value = 1)
x = np.array[r_i, 0, 0, ]

frm = ttk.Frame(root, padding=10)
frm.grid()

ttk.Label(frm, text="Initial altitude").grid(column=0, row=0)
a_i_Entry = ttk.Entry(frm, textvariable = a_i)
a_i_Entry.grid(column = 1, row = 0)
a_i.trace_add("write", altitude_changed)
ttk.Label(frm, text="km").grid(column = 2, row = 0)

ttk.Label(frm, text="Target altitude").grid(column=0, row=1)
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

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize = [9, 9], layout = "tight")
fig.set_facecolor("midnightblue")
fig.set_frameon(True)
fig.set_edgecolor("darkgrey")
ax.spines['polar'].set_visible(False)
ax.set_facecolor("midnightblue")
planet = Circle((0,0), R, transform=ax.transData._b, facecolor = "teal", edgecolor = "dimgrey")
initial_orbit = Circle((0,0), transform=ax.transData._b, edgecolor = "red", facecolor = '#00000000')
final_orbit = Circle((0,0), transform=ax.transData._b, edgecolor = "skyblue", facecolor = '#00000000')
transfer_orbit = Line2D([0,1],[0,1], color = "limegreen")
transfer_orbit_trace = Line2D([0,1],[0,1], color = "limegreen", linestyle = 'dashed')
ship, = ax.plot(0, 0, color = "white", marker = "^", markersize = 15)
station, = ax.plot(0, 0, color = "skyblue", marker = "D", markersize = 12)
rendezvous_point, = ax.plot(0, 0, color = "orange", marker = "X", markersize = 12)

ax.add_patch(planet)
ax.add_patch(initial_orbit)
ax.add_patch(final_orbit)
ax.add_line(transfer_orbit)
ax.add_line(transfer_orbit_trace)

ax.set_rticks([])			# No radial ticks
ax.set_rgrids([])			# No radial grids
ax.set_thetagrids([])		# No angular grids
ax.set_theta_offset(np.pi)	# put θ = 0 on the left
ax.set_theta_direction(-1)	# make positive θ clockwise
ax.grid(False)

canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.

Hohmann_transfer(r_i, r_f)

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
