#GROMACS Tutorials

Some [GROMACS](http://www.gromacs.org) tutorials for beginners. It's not
necessary to do the tutorials in order, but the first two
tutorials are essential before going on to the others, as the structure file
(`methane.pdb`) and topology file (`topol.top`) for methane from tutorial 2 are
used in all subsequent tutorials. The tutorials are designed for GROMACS version
5.1. If you are using an older version, some of the commands or parameters may
have changed. Note especially that the pull code for umbrella sampling has
changed since 5.0 and older releases.

I assume you have some working knowledge of the command line. Specifically, you
should know how to make directories, change into them, and edit text files.
Obviously I also assume you have GROMACS installed on a machine available to you.

Throughout the tutorials we'll be using OPLS methane and TIP4PEW water.

##Contents

1. [Water](https://github.com/wesbarnett/gromacs-tutorials/blob/master/1_tip4pew_water/README.md) - Basics of setting up a simulation. Find out the
   density of TIP4PEW water.
2. [One methane in water](https://github.com/wesbarnett/gromacs-tutorials/blob/master/2_methane_in_water/README.md) - How to create a topology file
   for a molecule and solvate it. Get the radial distribution function.
3. [Several methanes in water](https://github.com/wesbarnett/gromacs-tutorials/blob/master/3_methanes_in_water/README.md) - How to put multiple
   solutes into a system. Get the methane-methane potential of mean force.
4. [Free energy of solvation of methane](https://github.com/wesbarnett/gromacs-tutorials/blob/master/4_methane_fe/README.md) - How to do a free energy
   simulation when coupling a molecule. Use MBAR to get the result.
5. [Umbrella sampling](https://github.com/wesbarnett/gromacs-tutorials/blob/master/5_umbrella/README.md) - Get methane-methane PMF from umbrella sampling using pull
   code.

## Links

Some of the other software that I use in these tutorials that you may find
useful are:

* [alchemical-analysis](https://github.com/MobleyLab/alchemical-analysis)
* [Avogadro](http://avogadro.cc/wiki/Main_Page)
* [gnuplot](http://www.gnuplot.info/)
* [vmd](http://www.ks.uiuc.edu/Research/vmd/)

