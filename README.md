# FoN-Orbital-Calculator
An inter-orbit transfer and rendezvous tool intended for the game Flight of Nova, written in Python.

# Installation

## Windows/Mac
If you do not have Python 3 on your system, install it from https://www.python.org/downloads.

## Linux
Python is included in almost all distributions, but not all of these have version 3. To check, type `python --version` or `python3 --version` into the terminal, and if the command is not recognised, or the version is less than 3.12.x, install it according to your distribution or from the above link.

After completing the OS-specific installation of Python, some modules will need to be installed, preferably using `pip install -r requirements.txt`. Finally, download the repository and extract it.

# Usage
Run `FoN-orbit-calc.py`. Use the tabs at the top to select the type of maneuver to calculate.

## Hohmann Transfer Calculator
This section of the program assumes circular initial orbits for both your craft and your target. If not at a station, make sure your orbit is as circular as possible.\
Enter the altitude of your current orbit and the altitude of your target. The "Delta longitude" is the longitude difference between you and your target at which you should perform your burn. Positive values mean the target should be ahead of you in the orbit, negative values mean it should be behind. This allows the target to 'catch up' with you if you are coming from a higher, slower orbit, or for you to catch up with it if you are coming from a faster orbit. \
When the difference between your longitude and the station's longitude is at the right value (or if the distance between you and the station matches the value given by the calculator), burn so that your orbit just touches that of the target. So if you are ascending, you should match your apoapsis with the target's altitude; if you're descending, match your periapsis with target's altitude. Then, wait :)

### Glossary of parameters
Burn delta longitude: The difference in longitude between your ship and your target at the point when you should burn. Positive values mean the station should be 'ahead' of you in its orbit (its longitude should be higher than yours), while negative values mean it should be 'behind' (its longitude should be lower than yours). \
Distance to target: The distance between your ship and your target at the point when you should burn.
Transfer duration: The time it will take to reach your target, if you make the currently displayed burn.

## Custom Transfer Calculator
In this section, you can use the sliders to specify your burn direction and magnitude. It is then your responsibility to make sure your trajectory reaches the target altitude and does not collide with the planet. The horizontal slider controls the prograde/retrograde (forward/backward) component of the burn, and the vertical slider controls the radial out/radial in (towards/away from the planet) component. The ship image will give a rough idea of what orientation you should be in during your burn, where prograde is defined as to the right on the screen.

### Glossary of parameters
Burn AoE: The angle above horizontal in which you should point your ship during your burn. \
Burn deltaV: The amount of deltaV you will use, or equivalently, the amount by which your speed should have changed after your burn. \
Post burn FPA: The Flight Path Angle which you should aim to have immediately after your burn. \
Post burn speed: The speed you should try to have immediately after your burn.

## Reentry Calculator
This is an adapted version of the Custom Transfer Calculator with the "target" removed. Use this in the same way as the Custom Transfer Calculator to get a rough estimate of how far around the planet you will travel after your de-orbit burn. Atmospheric drag is not (yet!) factored into the calculation.
