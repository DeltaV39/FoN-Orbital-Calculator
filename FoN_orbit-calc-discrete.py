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

R = 6350    # Radius of NVA-3 in km
M = 5.9722e24
mu = 398600441800000 # Standard gravitational parameter of NVA-3 (assuming identical mass to earth)
scaled_G = scipy.constants.G/(1000**3)

def circular_orbit_dtheta(r):
    
    speed = np.sqrt((scaled_G*M)/np.power(r,3))
    
    return speed

# ~ def Hohmann_transfer(r_i, r_f):
    # ~ delta_longitude_degrees.set(180 * (1 - np.power((r_i + r_f)/(2*r_f), 3/2))) # the formula :O
    # ~ delta_longitude_degrees.set(round(delta_longitude_degrees.get() % (360*np.sign(delta_longitude_degrees.get())), 2)) # can be optimised
    # ~ if abs(delta_longitude_degrees.get()) > 180:
        # ~ delta_longitude_degrees.set(360 - abs(delta_longitude_degrees.get()))   # can DEFINITELY be optimised
    
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

def gravity_ODE(x):
    r = x[0,0]
    theta = x[0,1]
    dr = x[1,0]
    dtheta = x[1,1]
    d2r = r*np.power(dtheta,2) - (scaled_G*M)/np.power(r,2)
    d2theta = -(2*dr*dtheta)/r
    
    return np.array([[dr, dtheta], [d2r, d2theta]])

def velocity_verlet(ODE, s, timestep):
    v_mid = (s+ODE(s)*0.5*timestep)[1]
    x_end = s[0]+timestep*v_mid
    a_end = ODE(np.array([x_end, v_mid]))[1]
    v_end = v_mid + 0.5*timestep*a_end
    
    return(np.array([x_end, v_end], dtype = float))

def calc_trajectory():
    result = [initial_state]
    print(result)
    while (result[-1][0,1] < 2*np.pi) and (70+R < result[-1][0,0] < r_f):
        result.append(velocity_verlet(gravity_ODE, result[-1], 1))
    
    states_array = np.array(result)
    transposed = np.transpose(states_array, (1,2,0))
    # ~ print(transposed)
    return transposed[0]

def altitude_changed(*args):
    try:
        r_i = a_i.get() + R
        r_f = a_f.get() + R
        initial_state = np.array([[r_i, 0], [0, circular_orbit_dtheta(r_i)+tangential_deltav.get()]])
        print(initial_state)
    except TclError:
        print("You entered something strange")
        return
    
    # update plot with new altitudes
    trajectory = calc_trajectory()
    update_plot(trajectory)
    return
    
def burn_changed(*args):
    # update plot with new trajectory
    initial_state[1,0] = radial_deltav.get()
    initial_state[1,1] = circular_orbit_dtheta(r_i) + tangential_deltav.get()/r_i
    trajectory = calc_trajectory()
    update_plot(trajectory)
    # update ship sprite
    return

def update_plot(trajectory):
    initial_orbit.set(radius = r_i)
    final_orbit.set(radius = r_f)
    transfer_orbit.set_data(trajectory[1],trajectory[0])
    ship.set_data(0, r_i)
    # ~ station.set_data(np.radians(delta_longitude_degrees.get()), r_f)
    if abs(trajectory[0,-1] - r_f) < 100:
        rendezvous_point.set_data(trajectory[1,-1], r_f)
    ax.set_rmax(max(r_i, r_f)*1.1)
    
    canvas.draw()
    return

root = Tk()
root.wm_title("Inter-orbital Rendezvous Calculator")

a_i = DoubleVar(value = 185.0)
a_f = DoubleVar(value = 1058.0)
r_i = 185.0 + R
r_f = 1058.0 + R
delta_longitude_degrees = DoubleVar()
burn_distance = DoubleVar(value = 1000)
transfer_duration = DoubleVar(value = 1)
initial_state = np.array([[r_i, 0], [0, circular_orbit_dtheta(r_i)]]) # state vector (r, θ, dr, dθ)
print(initial_state)
radial_deltav = DoubleVar(value = 0)
tangential_deltav = DoubleVar(value = 0)

main_frame = ttk.Frame(root, padding=10)
main_frame.grid()

metrics_frame = ttk.Frame(main_frame, padding=10)
metrics_frame.grid(column = 0, row = 0)

ttk.Label(metrics_frame, text="Initial altitude").grid(column=0, row=0)
a_i_Entry = ttk.Entry(metrics_frame, textvariable = a_i)
a_i_Entry.grid(column = 1, row = 0)
a_i.trace_add("write", altitude_changed)
ttk.Label(metrics_frame, text="km").grid(column = 2, row = 0)

ttk.Label(metrics_frame, text="Target altitude").grid(column=0, row=1)
a_f_Entry = ttk.Entry(metrics_frame, textvariable = a_f)
a_f_Entry.grid(column = 1, row = 1)
a_f.trace_add("write", altitude_changed)
ttk.Label(metrics_frame, text="km").grid(column = 2, row = 1)

ttk.Label(metrics_frame, text = "Delta Longitude").grid(column = 0, row = 2)
delta_longitude_degrees_label = ttk.Label(metrics_frame, textvariable = delta_longitude_degrees)
delta_longitude_degrees_label.grid(column = 1, row = 2)
ttk.Label(metrics_frame, text="⁰").grid(column = 2, row = 2)

ttk.Label(metrics_frame, text = "Distance to target").grid(column = 0, row = 3)
burn_distance_label = ttk.Label(metrics_frame, textvariable = burn_distance)
burn_distance_label.grid(column = 1, row = 3)
ttk.Label(metrics_frame, text="km").grid(column = 2, row = 3)

ttk.Label(metrics_frame, text = "Transfer duration").grid(column = 0, row = 4)
transfer_duration_label = ttk.Label(metrics_frame, textvariable = transfer_duration)
transfer_duration_label.grid(column = 1, row = 4)
ttk.Label(metrics_frame, text="min").grid(column = 2, row = 4)

scales_frame = ttk.Frame(main_frame, padding=10)
scales_frame.grid(column = 0, row = 1)

radial_scale = ttk.Scale(
    scales_frame, 
    orient=VERTICAL, 
    length=200, 
    from_=5.0, 
    to=-5.0, 
    variable=radial_deltav,
    command=burn_changed)
radial_scale.grid(column = 0, row = 0)

tangential_scale = ttk.Scale(
    scales_frame, 
    orient=HORIZONTAL, 
    length=200, 
    from_=-5.0, 
    to=5.0, 
    variable=tangential_deltav,
    command=burn_changed)
tangential_scale.grid(column = 1, row = 1)

# ~ set up canvas for orbit plots
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize = [9, 9], layout = "tight")
fig.set_facecolor("midnightblue")
fig.set_frameon(True)
fig.set_edgecolor("darkgrey")
ax.spines['polar'].set_visible(False)
ax.set_facecolor("midnightblue")
theta = np.linspace(0, 2*np.pi, num = 1000)
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

ax.set_rticks([])           # No radial ticks
ax.set_rgrids([])           # No radial grids
ax.set_thetagrids([])       # No angular grids
ax.set_theta_offset(np.pi)  # put θ = 0 on the left
ax.set_theta_direction(-1)  # make positive θ clockwise
ax.grid(False)              # No grids

canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.

canvas.get_tk_widget().grid(column = 4, row = 0, padx = 4, pady = 4)

trajectory = calc_trajectory()
update_plot(trajectory)

def on_closing():
    root.destroy()
    quit()

root.protocol("WM_DELETE_WINDOW", on_closing)
root.mainloop()

# ~ a_i = float(input("Enter initial altitude in km: ")) *1000
# ~ a_f = float(input("Enter final altitude in km: ")) * 1000

# ~ delta_longitude_degrees = 180 * (1 - np.power((2*R + a_i + a_f)/(2*R+2*a_f), 3/2))

# ~ print(f"Burn at relative longitude of {delta_longitude_degrees}")
