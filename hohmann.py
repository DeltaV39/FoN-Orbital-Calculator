from tkinter import *
from tkinter import ttk

import numpy as np
import matplotlib.pyplot as plt
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure
from matplotlib.patches import Circle

class HohmannCalculator(ttk.Frame):

    def __init__(self, parent):
        super().__init__(parent)
        
        self.planet_R = 6350000 # Radius of NVA-3
        self.mu = 398600441800000 # Standard gravitational parameter of NVA-3 (assuming identical mass to earth)

        self.a_i = DoubleVar(value = 185.0)
        self.a_f = DoubleVar(value = 1058.0)
        self.delta_longitude_degrees = DoubleVar()
        self.burn_distance = DoubleVar(value = 1000)
        self.transfer_duration = DoubleVar(value = 1)

        self.frm = ttk.Frame(self, padding=10)
        self.frm.grid()

        ttk.Label(self.frm, text="Initial altitude").grid(column=0, row=0)
        self.a_i_Entry = ttk.Entry(self.frm, textvariable = self.a_i)
        self.a_i_Entry.grid(column = 1, row = 0)
        self.a_i.trace_add("write", self.altitude_changed)
        ttk.Label(self.frm, text="km").grid(column = 2, row = 0)

        ttk.Label(self.frm, text="Final altitude").grid(column=0, row=1)
        self.a_f_Entry = ttk.Entry(self.frm, textvariable = self.a_f)
        self.a_f_Entry.grid(column = 1, row = 1)
        self.a_f.trace_add("write", self.altitude_changed)
        ttk.Label(self.frm, text="km").grid(column = 2, row = 1)

        ttk.Label(self.frm, text = "Delta Longitude").grid(column = 0, row = 2)
        self.delta_longitude_degrees_label = ttk.Label(self.frm, textvariable = self.delta_longitude_degrees)
        self.delta_longitude_degrees_label.grid(column = 1, row = 2)
        ttk.Label(self.frm, text="⁰").grid(column = 2, row = 2)

        ttk.Label(self.frm, text = "Distance to target").grid(column = 0, row = 3)
        self.burn_distance_label = ttk.Label(self.frm, textvariable = self.burn_distance)
        self.burn_distance_label.grid(column = 1, row = 3)
        ttk.Label(self.frm, text="km").grid(column = 2, row = 3)

        ttk.Label(self.frm, text = "Transfer duration").grid(column = 0, row = 4)
        self.transfer_duration_label = ttk.Label(self.frm, textvariable = self.transfer_duration)
        self.transfer_duration_label.grid(column = 1, row = 4)
        ttk.Label(self.frm, text="min").grid(column = 2, row = 4)

        # set up self.canvas for orbit plots

        self.theta = np.linspace(0, 2*np.pi, num = 1000)
        self.r_i = self.a_i.get()*1000 + self.planet_R
        self.r_f = self.a_f.get()*1000 + self.planet_R

        # elliptical orbit parameters

        self.orbit_param_c = -(self.r_f-self.r_i)/2

        self.orbit_param_a = (self.r_i+self.r_f)/2

        self.r_transfer = (-self.orbit_param_a**2 +self.orbit_param_c**2)/(-self.orbit_param_a + self.orbit_param_c*np.cos(self.theta))

        self.fig, self.ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize = [9, 9], layout = "tight")
        self.fig.set_facecolor("midnightblue")
        self.fig.set_frameon(True)
        self.fig.set_edgecolor("darkgrey")
        self.ax.spines['polar'].set_visible(False)
        self.ax.set_facecolor("midnightblue")
        self.planet = Circle((0,0), self.planet_R, transform=self.ax.transData._b, facecolor = "teal", edgecolor = "dimgrey")
        self.ax.add_artist(self.planet)
        self.initial_orbit, = self.ax.plot(self.theta, np.full_like(self.theta, self.r_i), color = "red")
        self.final_orbit, = self.ax.plot(self.theta, np.full_like(self.theta, self.r_f), color = "skyblue")
        self.transfer_orbit, = self.ax.plot(self.theta[:self.theta.size//2], self.r_transfer[:self.theta.size//2], color = "limegreen")
        self.transfer_orbit_trace, = self.ax.plot(self.theta[self.theta.size//2:], self.r_transfer[self.theta.size//2:], color = "limegreen", linestyle = 'dashed')
        self.ship, = self.ax.plot(0, self.r_i, color = "white", marker = "^", markersize = 15)
        self.station, = self.ax.plot(np.radians(self.delta_longitude_degrees.get()), self.r_f, color = "skyblue", marker = "D", markersize = 12)
        self.rendezvous_point, = self.ax.plot(np.pi, self.r_f, color = "orange", marker = "X", markersize = 12)

        self.ax.set_rmax(max(self.r_i, self.r_f)*1.1)
        self.ax.set_rticks([])   # No radial ticks
        self.ax.set_rgrids([])   # No radial grids
        self.ax.set_theta_offset(np.pi)
        self.ax.set_theta_direction(-1)
        self.ax.set_thetagrids([])
        self.ax.grid(False)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)  # A tk.DrawingArea.

        self.recalc_burn_params()

        self.canvas.get_tk_widget().grid(column = 4, row = 0, padx = 4, pady = 4)

    def ellipse_TA_from_radius(self, r, ecc):
        return(np.acos((self.orbit_param_a*(1-ecc**2))/(r*ecc)))

    def recalc_burn_params(self):
        self.delta_longitude_degrees.set(180 * (1 - np.power((self.r_i + self.r_f)/(2*self.r_f), 3/2))) # the formula :O
        self.delta_longitude_degrees.set(round(self.delta_longitude_degrees.get() % (360*np.sign(self.delta_longitude_degrees.get())), 2)) # can be optimised
        if abs(self.delta_longitude_degrees.get()) > 180:
            self.delta_longitude_degrees.set(360 - abs(self.delta_longitude_degrees.get()))   # can DEFINITELY be optimised
        
        self.burn_distance.set(round(np.sqrt(self.r_i**2+self.r_f**2 - 2*self.r_i*self.r_f*np.cos(np.radians(self.delta_longitude_degrees.get())))/1000, 2))
        
        self.transfer_duration.set(round((np.pi*np.sqrt((((self.r_i+self.r_f)/2)**3)/self.mu))/60, 1))
        
        self.orbit_param_c = -(self.r_f-self.r_i)/2
        
        self.orbit_param_a = (self.r_i+self.r_f)/2
        
        self.update_plot()

    def update_plot(self):
        self.theta = np.linspace(0, 2*np.pi, num = 1000)
        self.r_transfer = (-self.orbit_param_a**2 +self.orbit_param_c**2)/(-self.orbit_param_a + self.orbit_param_c*np.cos(self.theta))
        
        self.initial_orbit.set_data(self.theta, np.full_like(self.theta, self.r_i))
        self.final_orbit.set_data(self.theta, np.full_like(self.theta, self.r_f))
        self.transfer_orbit.set_data(self.theta[:self.theta.size//2], self.r_transfer[:self.theta.size//2])
        self.transfer_orbit_trace.set_data(self.theta[self.theta.size//2:], self.r_transfer[self.theta.size//2:])
        self.ship.set_data([0], [self.r_i])
        self.station.set_data([np.radians(self.delta_longitude_degrees.get())], [self.r_f])
        self.rendezvous_point.set_data([np.pi], [self.r_f])
        self.ax.set_rmax(max(self.r_i, self.r_f)*1.1)
        
        self.canvas.draw()
        return

    def altitude_changed(self, *args):
        try:
            self.r_i = self.a_i.get()*1000 + self.planet_R
            self.r_f = self.a_f.get()*1000 + self.planet_R
        
        except TclError:
            return
        
        self.recalc_burn_params()
        return
