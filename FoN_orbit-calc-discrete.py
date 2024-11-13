import math
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
M = 5.9722e24   # Mass of NVA-3 in kg
mu = 398600441800000 # Standard gravitational parameter of NVA-3 (assuming identical mass to earth)
scaled_G = scipy.constants.G/(1000**3)  # Scaled to work with km 

def circular_orbit_speed(r):
    
    return np.sqrt((scaled_G*M)/np.power(r,3))

def circular_orbit_progress_angle(r, t):
    
    return math.tau*(t/math.sqrt((4*(math.pi**2)/(scaled_G*M))*r**3))

def gravity_ODE(x):
    r = x[0,0]
    theta = x[0,1]
    dr = x[1,0]
    dtheta = x[1,1]
    d2r = r*(dtheta**2) - (scaled_G*M)/(r**2)
    d2theta = -(2*dr*dtheta)/r
    
    return np.array([[dr, dtheta], [d2r, d2theta]])

def velocity_verlet(ODE, s, timestep):
    v_mid = (s+ODE(s)*0.5*timestep)[1]
    x_end = s[0]+timestep*v_mid
    a_end = ODE(np.array([x_end, v_mid]))[1]
    v_end = v_mid + 0.5*timestep*a_end
    
    return(np.array([x_end, v_end], dtype = float))

def calc_trajectory():
    global encounter_exists, initial_state
    result = [initial_state]
    while(
        (result[-1][0,1] < math.tau)
        and (result[-1][0,0] > R)
        and (abs(result[-1][0,0] - r_f) > 10)
        and len(result) < 7200
        ):
        result.append(velocity_verlet(gravity_ODE, result[-1], 1))
    
    if abs(result[-1][0,0] - r_f) < 10:
        encounter_exists = True
    else:
        encounter_exists = False
    states_array = np.array(result)
    transposed = np.transpose(states_array, (1,2,0))
    return transposed[0]

def altitude_changed(*args):
    global r_i, r_f, initial_state
    try:
        r_i = a_i.get() + R
        r_f = a_f.get() + R
        initial_state = np.array([[r_i, 0], [0, circular_orbit_speed(r_i)+tangential_deltav.get()]])
    except TclError:
        print("You entered something strange")
        return
    
    # update plot with new altitude and trajectory
    trajectory =  calc_trajectory()
    update_parameters(trajectory)
    update_plot(trajectory)
    return
    
def burn_changed(*args):
    # update plot with new trajectory
    initial_state[1,0] = radial_deltav.get()
    initial_state[1,1] = circular_orbit_speed(r_i) + tangential_deltav.get()/r_i
    trajectory =  calc_trajectory()
    update_parameters(trajectory)
    update_plot(trajectory)
    # update ship sprite
    return

def update_parameters(trajectory):
    global delta_longitude_radians
    if encounter_exists:
        transfer_duration_secs = trajectory.shape[1]    # this works because each step is 1 second after the last
        transfer_duration_mins.set(round(transfer_duration_secs/60,1)) 
        # set delta longitude degrees and burn distance
        delta_longitude_radians = trajectory[1,-1]-circular_orbit_progress_angle(r_f, transfer_duration_secs)
        delta_longitude_degrees.set(round(math.degrees(delta_longitude_radians), 1))
        
    else:
        delta_longitude_degrees.set("N/A")
        burn_distance.set("N/A")
        transfer_duration_mins.set("∞")
    try:
        burn_AoE.set(round(math.degrees(math.atan(radial_deltav.get()/tangential_deltav.get())),1))
        burn_AoE.set(math.copysign(burn_AoE.get(), radial_deltav.get()))
    except ZeroDivisionError:
        burn_AoE.set(90)
        burn_AoE.set(math.copysign(burn_AoE.get(), radial_deltav.get()))
    burn_deltav.set(round(np.sqrt(radial_deltav.get()**2 + tangential_deltav.get()**2)*1000,1))
    burn_distance.set(round((np.sqrt(r_i**2+r_f**2 - 2*r_i*r_f*math.cos(delta_longitude_radians))), 1))
    try:
        after_FPA.set(round(math.degrees(math.atan(radial_deltav.get()/(circular_orbit_speed(r_i)*r_i+tangential_deltav.get()))),1))
    except ZeroDivisionError:
        after_FPA.set(90)
    
    after_speed.set(round((math.sqrt((radial_deltav.get())**2+(circular_orbit_speed(r_i)*r_i+tangential_deltav.get())**2))*1000,1))
    
    # set appropriate ship sprite for burn AoE
    if tangential_deltav.get() >= 0:
        sprite_num = int(round(-burn_AoE.get()/15, 0))+6
    else:
        sprite_num = int(round(burn_AoE.get()/15, 0))+18
    
    ship_sprite_label['image'] = sprites[sprite_num]
    return

def update_plot(trajectory):
    initial_orbit.set(radius = r_i)
    final_orbit.set(radius = r_f)
    transfer_orbit.set_data(trajectory[1],trajectory[0])
    ship.set_data(0, r_i)
    if encounter_exists:
        rendezvous_point.set_data(trajectory[1,-1], r_f)    # set orange cross at the end of the trajectory
        rendezvous_point.set_visible(True)
        # set station location
        station.set_data(delta_longitude_radians, r_f)   # put station where it would be at the time of the burn
        station.set_visible(True)
    else:
        rendezvous_point.set_visible(False)
        station.set_visible(False)
    ax.set_rmax(max(r_i, r_f)*1.1)
    
    canvas.draw()
    return

root = Tk()
root.wm_title("Inter-orbital Rendezvous Calculator")

a_i = DoubleVar(value = 185.0)          # intial altitude
a_f = DoubleVar(value = 1058.0)         # final (target) altitude
radial_deltav = DoubleVar(value = 0)    # radial out component of burn
tangential_deltav = DoubleVar(value = 0)# prograde component of burn
delta_longitude_degrees = DoubleVar()   # longitude difference to target during burn
delta_longitude_radians = 0             # radians version for calculations
burn_distance = DoubleVar(value = 1000) # distance to target during burn 
burn_AoE = DoubleVar(value = 0)         # AoE while burning
burn_deltav = DoubleVar(value = 0)      # deltaV used during the burn
after_FPA = DoubleVar(value = 0)        # post-burn FPA
after_speed = DoubleVar(value = 8000)   # post-burn speed

transfer_duration_mins = DoubleVar(value = 1)
encounter_exists = False

r_i = 185.0 + R
r_f = 1058.0 + R
initial_state = np.array([[r_i, 0], [0, circular_orbit_speed(r_i)]]) # state vector (r, θ, dr, dθ)

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

ttk.Label(metrics_frame, text = "Burn delta longitude").grid(column = 0, row = 2)
delta_longitude_degrees_label = ttk.Label(metrics_frame, textvariable = delta_longitude_degrees)
delta_longitude_degrees_label.grid(column = 1, row = 2)
ttk.Label(metrics_frame, text="⁰").grid(column = 2, row = 2)

ttk.Label(metrics_frame, text = "Distance to target").grid(column = 0, row = 3)
burn_distance_label = ttk.Label(metrics_frame, textvariable = burn_distance)
burn_distance_label.grid(column = 1, row = 3)
ttk.Label(metrics_frame, text="km").grid(column = 2, row = 3)

ttk.Label(metrics_frame, text = "Burn AoE").grid(column = 0, row = 4)
burn_distance_label = ttk.Label(metrics_frame, textvariable = burn_AoE)
burn_distance_label.grid(column = 1, row = 4)
ttk.Label(metrics_frame, text="⁰").grid(column = 2, row = 4)

ttk.Label(metrics_frame, text = "Burn deltaV").grid(column = 0, row = 5)
burn_distance_label = ttk.Label(metrics_frame, textvariable = burn_deltav)
burn_distance_label.grid(column = 1, row = 5)
ttk.Label(metrics_frame, text="m/s").grid(column = 2, row = 5)

ttk.Label(metrics_frame, text = "Post burn FPA").grid(column = 0, row = 6)
burn_distance_label = ttk.Label(metrics_frame, textvariable = after_FPA)
burn_distance_label.grid(column = 1, row = 6)
ttk.Label(metrics_frame, text="⁰").grid(column = 2, row = 6)

ttk.Label(metrics_frame, text = "Post burn speed").grid(column = 0, row = 7)
burn_distance_label = ttk.Label(metrics_frame, textvariable = after_speed)
burn_distance_label.grid(column = 1, row = 7)
ttk.Label(metrics_frame, text="m/s").grid(column = 2, row = 7)

ttk.Label(metrics_frame, text = "Transfer duration").grid(column = 0, row = 8)
transfer_duration_mins_label = ttk.Label(metrics_frame, textvariable = transfer_duration_mins)
transfer_duration_mins_label.grid(column = 1, row = 8)
ttk.Label(metrics_frame, text="min").grid(column = 2, row = 8)

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

ship_sprite_label = ttk.Label(scales_frame, justify = 'center')
ship_sprite_label.grid(column = 1, row = 0)

sprites = []
# import right-facing ship sprites
for i in range(13):
    sprites.append(PhotoImage(file=f"ship_sprites/ship_right{i}.png"))

# import left-facing ship sprites
for i in range(11):
    sprites.append(PhotoImage(file=f"ship_sprites/ship_left{i+1}.png"))

# now sprites[x] represents the xth clockwise rotation of the ship by 15
# degrees, starting vertically up at 0.
ship_sprite_label['image'] = sprites[6]

# ~ set up canvas for orbit plots
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize = [9, 9], layout = "tight")
fig.set_facecolor("midnightblue")
fig.set_frameon(True)
fig.set_edgecolor("darkgrey")
ax.spines['polar'].set_visible(False)
ax.set_facecolor("midnightblue")
theta = np.linspace(0, math.tau, num = 1000)
planet = Circle((0,0), R, transform=ax.transData._b, facecolor = "teal", edgecolor = "dimgrey")
initial_orbit = Circle((0,0), transform=ax.transData._b, edgecolor = "red", facecolor = '#00000000')
final_orbit = Circle((0,0), transform=ax.transData._b, edgecolor = "skyblue", facecolor = '#00000000')
transfer_orbit = Line2D([0,1],[0,1], color = "limegreen")
ship, = ax.plot(0, 0, color = "white", marker = "^", markersize = 15)
station, = ax.plot(0, 0, color = "skyblue", marker = "D", markersize = 12)
rendezvous_point, = ax.plot(0, 0, color = "orange", marker = "X", markersize = 12)

# make markers draw above orbit lines
ship.set_zorder(9)
station.set_zorder(9)
rendezvous_point.set_zorder(9)

ax.add_patch(planet)
ax.add_patch(initial_orbit)
ax.add_patch(final_orbit)
ax.add_line(transfer_orbit)

ax.set_rticks([])           # No radial ticks
ax.set_rgrids([])           # No radial grids
ax.set_thetagrids([])       # No angular grids
ax.set_theta_offset(math.pi)  # put θ = 0 on the left
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
