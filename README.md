# KSP deltaV plots

Quick reference for stage design in KSP. Refer to the figures in the
dvplots2 folder -- the dvplots folder is a shittier version.

## How to read the figure

Figures are titled `engine-[atm/vac].png`. `atm/vac` refers to the engine's
performance at Kerbin surface or in a vacuum.

Each figure encodes possible ship configurations as a function of the cargo
mass (`topm`, the x-axis), and the desired delta-v of the entire stage (`dv`,
the y-axis).

There are three sets of lines on each figure:

1. In blue (with legend at right) are the possible ship configurations in terms
of possible number of fuel tanks per engine. Darker lines include more fuel
tanks. Use the legend title to figure out which fuel tanks are being modeled.
Numbers above the lines indicate the ratio of cargo to stage mass at that
configuration.

2. In light red/purple, are lines of constant cargo/stage mass.

3. In yellow/green gradient (with legend at top) are the thrust/weight ratios
for each configuration (referring either to Kerbin or Mun gravity).

