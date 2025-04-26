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

class ReentryCalculator(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.planet_R = 6350    # Radius of NVA-3 in km
        self.planet_M = 5.9722e24   # Mass of NVA-3 in kg
        self.mu = 398600441800000 # Standard gravitational parameter of NVA-3 (assuming identical mass to earth)
        self.scaled_G = scipy.constants.G/(1000**3)  # Scaled to work with km 

        self.a_i = DoubleVar(value = 185.0)          # intial altitude
        self.radial_deltav = DoubleVar(value = 0)    # radial out component of burn
        self.tangential_deltav = DoubleVar(value = 0)# prograde component of burn
        self.burn_AoE = DoubleVar(value = 0)         # AoE while burning
        self.burn_deltav = DoubleVar(value = 0)      # deltaV used during the burn
        self.after_FPA = DoubleVar(value = 0)        # post-burn FPA
        self.after_speed = DoubleVar(value = 8000)   # post-burn speed

        self.transfer_duration_mins = DoubleVar(value = 1)
        self.planetfall_distance = DoubleVar(value = 1)

        self.r_i = 185.0 + self.planet_R
        self.r_f = 185.0 + self.planet_R
        self.initial_state = np.array([[self.r_i, 0], [0, self.circular_orbit_speed(self.r_i)]]) # state vector [[r, θ], [dr, dθ]]

        self.main_frame = ttk.Frame(self, padding=10)
        self.main_frame.grid()

        self.metrics_frame = ttk.Frame(self.main_frame, padding=10)
        self.metrics_frame.grid(column = 0, row = 0)
        
        self.invalid_outline = ttk.Style()
        self.invalid_outline.configure("invalid.TEntry", fieldbackground='IndianRed1')
        ttk.Label(self.metrics_frame, text="Initial altitude").grid(column=0, row=0)
        self.a_i_Entry = ttk.Entry(self.metrics_frame, textvariable = self.a_i)
        self.a_i_Entry.grid(column = 1, row = 0)
        self.a_i.trace_add("write", self.altitude_changed)
        ttk.Label(self.metrics_frame, text="km").grid(column = 2, row = 0)

        ttk.Label(self.metrics_frame, text = "Burn AoE").grid(column = 0, row = 4)
        self.burn_AOE_label = ttk.Label(self.metrics_frame, textvariable = self.burn_AoE)
        self.burn_AOE_label.grid(column = 1, row = 4)
        ttk.Label(self.metrics_frame, text="⁰").grid(column = 2, row = 4)

        ttk.Label(self.metrics_frame, text = "Burn deltaV").grid(column = 0, row = 5)
        self.burn_deltaV_label = ttk.Label(self.metrics_frame, textvariable = self.burn_deltav)
        self.burn_deltaV_label.grid(column = 1, row = 5)
        ttk.Label(self.metrics_frame, text="m/s").grid(column = 2, row = 5)

        ttk.Label(self.metrics_frame, text = "Post burn FPA").grid(column = 0, row = 6)
        self.burn_FPA_label = ttk.Label(self.metrics_frame, textvariable = self.after_FPA)
        self.burn_FPA_label.grid(column = 1, row = 6)
        ttk.Label(self.metrics_frame, text="⁰").grid(column = 2, row = 6)

        ttk.Label(self.metrics_frame, text = "Post burn speed").grid(column = 0, row = 7)
        self.burn_speed_label = ttk.Label(self.metrics_frame, textvariable = self.after_speed)
        self.burn_speed_label.grid(column = 1, row = 7)
        ttk.Label(self.metrics_frame, text="m/s").grid(column = 2, row = 7)

        ttk.Label(self.metrics_frame, text = "Time to planetfall").grid(column = 0, row = 8)
        self.transfer_duration_mins_label = ttk.Label(self.metrics_frame, textvariable = self.transfer_duration_mins)
        self.transfer_duration_mins_label.grid(column = 1, row = 8)
        ttk.Label(self.metrics_frame, text="min").grid(column = 2, row = 8)

        ttk.Label(self.metrics_frame, text = "Distance to planetfall").grid(column = 0, row = 9)
        self.planetfall_distance_label = ttk.Label(self.metrics_frame, textvariable = self.planetfall_distance)
        self.planetfall_distance_label.grid(column = 1, row = 9)
        ttk.Label(self.metrics_frame, text="km").grid(column = 2, row = 9)

        self.scales_frame = ttk.Frame(self.main_frame, padding=10)
        self.scales_frame.grid(column = 0, row = 1)

        self.radial_scale = ttk.Scale(
            self.scales_frame,
            orient=VERTICAL,
            length=200,
            from_=5.0,
            to=-5.0,
            variable=self.radial_deltav,
            command=self.burn_changed)
        self.radial_scale.grid(column = 0, row = 0)

        self.tangential_scale = ttk.Scale(
            self.scales_frame,
            orient=HORIZONTAL,
            length=200,
            from_=-5.0,
            to=5.0,
            variable=self.tangential_deltav,
            command=self.burn_changed)
        self.tangential_scale.grid(column = 1, row = 1)

        self.ship_sprite_label = ttk.Label(self.scales_frame, justify = 'center')
        self.ship_sprite_label.grid(column = 1, row = 0)

        self.sprites = []
        # import right-facing ship sprites
        for i in range(13):
            self.sprites.append(PhotoImage(file=f"ship_sprites/ship_right{i}.png"))

        # import left-facing ship sprites
        for i in range(11):
            self.sprites.append(PhotoImage(file=f"ship_sprites/ship_left{i+1}.png"))

        # now sprites[x] represents the xth clockwise rotation of the ship by 15
        # degrees, starting vertically up at 0.
        self.ship_sprite_label['image'] = self.sprites[6]

        # set up canvas for orbit plots
        self.fig, self.ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize = [9, 9], layout = "tight")
        self.fig.set_facecolor("midnightblue")
        self.fig.set_frameon(True)
        self.fig.set_edgecolor("darkgrey")
        self.ax.spines['polar'].set_visible(False)
        self.ax.set_facecolor("midnightblue")
        self.theta = np.linspace(0, math.tau, num = 1000)
        self.planet = Circle((0,0), self.planet_R, transform=self.ax.transData._b, facecolor = "teal", edgecolor = "dimgrey")
        self.initial_orbit = Circle((0,0), transform=self.ax.transData._b, edgecolor = "red", facecolor = '#00000000')
        self.transfer_orbit = Line2D([0,1],[0,1], color = "limegreen")
        self.ship, = self.ax.plot(0, 0, color = "white", marker = "^", markersize = 15)
        self.planetfall_point, = self.ax.plot(0, 0, color = "lime", marker = "X", markersize = 12)

        # make markers draw above orbit lines
        self.ship.set_zorder(9)
        self.planetfall_point.set_zorder(7)

        self.ax.add_patch(self.planet)
        self.ax.add_patch(self.initial_orbit)
        self.ax.add_line(self.transfer_orbit)

        self.ax.set_rticks([])               # No radial ticks
        self.ax.set_rgrids([])               # No radial grids
        self.ax.set_thetagrids([])           # No angular grids
        self.ax.set_theta_offset(math.pi)    # put θ = 0 on the left
        self.ax.set_theta_direction(-1)      # make positive θ clockwise
        self.ax.grid(False)                  # No grids

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)  # A tk.DrawingArea.

        self.canvas.get_tk_widget().grid(column = 4, row = 0, padx = 4, pady = 4)

        self.trajectory = self.calc_trajectory()
        self.update_plot()

    def circular_orbit_speed(self, r):
        return np.sqrt((self.scaled_G*self.planet_M)/np.power(r,3))

    def gravity_ODE(self, x):
        r = x[0,0]
        self.theta = x[0,1]
        dr = x[1,0]
        dtheta = x[1,1]
        d2r = r*(dtheta**2) - (self.scaled_G*self.planet_M)/(r**2)
        d2theta = -(2*dr*dtheta)/r
        
        return np.array([[dr, dtheta], [d2r, d2theta]])

    def velocity_verlet(self, ODE, s, timestep):
        v_mid = (s+ODE(s)*0.5*timestep)[1]
        x_end = s[0]+timestep*v_mid
        a_end = ODE(np.array([x_end, v_mid]))[1]
        v_end = v_mid + 0.5*timestep*a_end
        
        return(np.array([x_end, v_end], dtype = float))

    # TODO: Add an actual array of times and adaptive timestep
    # TODO: Integrate the burn itself, otherwise integration is pointless
    def calc_trajectory(self):
        result = [self.initial_state]
        while(
            (result[-1][0,1] < math.tau)                # Don't integrate past a full orbit
            and (result[-1][0,0] > self.planet_R)       # Don't integrate past hitting the ground
            and len(result) < 7200                      # who would want a transfer to take over 2 hours lol
            ):
            result.append(self.velocity_verlet(self.gravity_ODE, result[-1], 1))
        
        states_array = np.array(result)
        transposed = np.transpose(states_array, (1,2,0))
        return transposed[0]    # Array of shape (2, 2, number_of_iterations), [[r_1 ... r_n, dr_1, ... dr_n], [theta_1, ... theta_n, dtheta1, ... dtheta_n]]

    def altitude_changed(self, *args):
        try:
            self.r_i = self.a_i.get() + self.planet_R
        except TclError:
            self.a_i_Entry.winfo_class()
            self.a_i_Entry.configure(style="invalid.TEntry")
            return
        self.a_i_Entry.configure(style="TEntry")
        self.initial_state = np.array([[self.r_i, 0], [0, self.circular_orbit_speed(self.r_i)+self.tangential_deltav.get()]])
        # update plot with new altitude and trajectory
        self.radial_deltav.set(0)
        self.tangential_deltav.set(0)
        self.trajectory =  self.calc_trajectory()
        self.update_parameters()
        self.update_plot()
        return
        
    def burn_changed(self, *args):
        # update plot with new trajectory
        self.initial_state[1,0] = self.radial_deltav.get()
        self.initial_state[1,1] = self.circular_orbit_speed(self.r_i) + self.tangential_deltav.get()/self.r_i
        self.trajectory =  self.calc_trajectory()
        self.update_parameters()
        self.update_plot()
        # update ship sprite
        return

    def update_parameters(self):
        transfer_duration_secs = self.trajectory.shape[1]    # this works because each step is 1 second after the last
        self.transfer_duration_mins.set(round(transfer_duration_secs/60,1))
        final_position = self.trajectory[:, -1]
        if final_position[0] <= self.planet_R:
            self.planetfall_distance.set(round(math.sqrt(self.r_i**2 + final_position[0]**2 - 2*self.r_i*final_position[0]*math.cos(final_position[1])), 1))
        else:
            self.planetfall_distance.set(math.inf)
        try:
            self.burn_AoE.set(round(math.degrees(math.atan(self.radial_deltav.get()/self.tangential_deltav.get())),1))
            self.burn_AoE.set(math.copysign(self.burn_AoE.get(), self.radial_deltav.get()))
        except ZeroDivisionError:
            self.burn_AoE.set(90)
            self.burn_AoE.set(math.copysign(self.burn_AoE.get(), self.radial_deltav.get()))
        self.burn_deltav.set(round(np.sqrt(self.radial_deltav.get()**2 + self.tangential_deltav.get()**2)*1000,1))
        try:
            self.after_FPA.set(round(math.degrees(math.atan(self.radial_deltav.get()/(self.circular_orbit_speed(self.r_i)*self.r_i+self.tangential_deltav.get()))),1))
        except ZeroDivisionError:
            self.after_FPA.set(90)
        
        self.after_speed.set(round((math.sqrt((self.radial_deltav.get())**2+(self.circular_orbit_speed(self.r_i)*self.r_i+self.tangential_deltav.get())**2))*1000,1))
        
        # set appropriate ship sprite for burn AoE
        if self.tangential_deltav.get() >= 0:
            sprite_num = int(round(-self.burn_AoE.get()/15, 0))+6
        else:
            sprite_num = (int(round(self.burn_AoE.get()/15, 0))+18) % 24
        
        self.ship_sprite_label['image'] = self.sprites[sprite_num]
        return

    def update_plot(self):
        self.initial_orbit.set(radius = self.r_i)
        self.transfer_orbit.set_data(self.trajectory[1],self.trajectory[0])
        self.ship.set_data([0], [self.r_i])
        if self.planetfall_distance.get() > 14000:
            self.planetfall_point.set_visible(False)
        else:
            self.planetfall_point.set_data([self.trajectory[1, -1]], [self.trajectory[0, -1]])
            self.planetfall_point.set_visible(True)
        self.ax.set_rmax(max(self.r_i, self.r_f, self.trajectory[0].max())*1.1)
        
        self.canvas.draw()
        return
