# FoN-Orbital-Calculator
An inter-orbit transfer and rendezvous tool intended for the game Flight of Nova, written in Python.

# Installation

## Windows/Mac
If you do not have Python 3 on your system, install it from https://www.python.org/downloads.

## Linux
Python is included in almost all distributions, but not all of these have version 3. To check, type `python3` into the terminal, and if the command is not recognised, install it according to your distribution or from the above link.

If Python 3 comes with your distribution, chances are it doesn't include the extra modules. These will need to be installed, either with pip or with your distribution's package management tool. In this case the modules are `python3-numpy`, `python3-tk`, and `python3-matplotlib`.

After completing the OS-specific installation of Python, download and run `FoN-orbit-calc.py`.

# Usage
This program assumes circular initial orbits for both your craft and your target.
Enter the altitude of your current orbit and the altitude of your target. The "Delta longitude" is the longitude difference between you and your target at which you should perform your burn. Positive values mean the target should be ahead of you in the orbit, negative values mean it should be behind. This allows the target to 'catch up' with you if you are coming from a higher, slower orbit. The reverse applies for a lower, faster orbit.
When the difference between your longitude and the station's longitude is at the right value, burn so that your orbit just touches that of the target. So if you are ascending, you should match your apoapsis with the target's altitude, and if you're descending, match your periapsis with target's altitude. Then, wait :)
