#GROMACS Tutorials

Some [GROMACS](http://www.gromacs.org) tutorials for beginners. Each tutorial
builds off of the previous, so it is recommended to do them in order.

I assume you have some working knowledge of the command line. Specifically, you
should know how to make directories, change into them, and edit text files.
Obviously I also assume you have GROMACS installed on a machine available to you.

Throughout the tutorials we'll be using OPLS methane and TIP4PEW water.

##Contents

1. [Water](1_tip4pew_water) - Basics of setting up a simulation. Find out the
   density of TIP4PEW water.
2. [One methane in water](2_methane_in_water) - How to create a topology file
   for a molecule and solvate it. Get the radial distribution function.
3. [Several methanes in water](3_methanes_in_water) - How to put multiple
   solutes into a system. Get the methane-methane potential of mean force.
4. [Free energy of solvation of methane](4_methane_fe) - How to do a free energy
   simulation when coupling a molecule. Use MBAR to get the result.
5. Umbrella sampling - Get methane-methane PMF from umbrella sampling using pull
   code. Coming soon.

## Links
Some of the other software that I use in these tutorials that you may find
useful are:

* [Avogadro](http://avogadro.cc/wiki/Main_Page)
* [gnuplot](http://www.gnuplot.info/)
* [vmd](http://www.ks.uiuc.edu/Research/vmd/)
