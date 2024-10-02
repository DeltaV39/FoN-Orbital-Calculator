import numpy as np
import tkinter as tk
import matplotlib as plt

R = 6350000	# Radius of NVA-3

a_i = float(input("Enter initial altitude in km: ")) *1000
a_f = float(input("Enter final altitude in km: ")) * 1000

delta_longitude = 180 * (1 - np.power((2*R + a_i + a_f)/(2*R+2*a_f), 3/2))

print(f"Burn at relative longitude of {delta_longitude}")
