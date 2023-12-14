# MotionAsym_Eye

Motion conditions: {cardinal, oblique}
Tilt angles: {± 0.5, ±1. ±2, ±4, ±8}

finalprocess14.m : script that counts all the micro saccades per motion condition and tilt angle. Plots the main comparison figures for VSS poster 2023. The micro saccades that are "counted" are only for trials with no blink and filtered for position within the "outside criteria". This uses a padding time of 0 (stimulus onset) to 300 (ms after stimulus offset).