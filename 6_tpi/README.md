[Return to main page](https://wesbarnett.github.io/gromacs-tutorials/)

Tutorial 6: Methane Chemical Potential from Test Particle Insertion
===================================================================

In this tutorial we'll randomly insert a united-atom methane into a box of water
thousands of times. Each time the particle is inserted, we can get its
interaction energy with the rest of the system. This can be used to get the
excess chemical potential of the methane from [Widom's
method](https://en.wikipedia.org/wiki/Widom_insertion_method).

Setup
-----

### Structure file

Copy `prd.gro` from the tutorial 1 where we ran a box of water and name it
`conf.gro`. If you haven't done that tutorial, you should do it now. You must
have a trajectory (.xtc) of water without the methane. In the file increase the
number of atoms (second line) by 1. Then add a single atom at the bottom of the
file before the box dimensions with residue `1CH4` and atom name `C`. You can
put whatever coordinates in, but the origin would be a good choice.

### Topology file

Reuse the topology file from tutorial 3. We're going to have to modify the file
a little bit so that we use a united atom methane instead of the all-atom
representation. Delete all of the hydrogens under `[ atoms ]`. For the carbon
change the atom type to `opls_066`. Change the charge to `0.000`. Change the
mass to . At the bottom of the file under `[ molecules ]` make sure the methane
is last.

### Parameter file

We're going to reuse the production file from tutorial 1 with some
modifications. Most of the parameters will be ignored by GROMACS, since this is
not a normal simulation. Set `md` to `tpi` and `nsteps` to `2000`. This tells
GROMACS to insert the test particle 2,000 times at each frame. Additionally
ensure that the reference temperature is the same as the original water
simulation (if you copied the file, it should be). Save this file in a folder
named `mdp` as `prd.mdp`.

Simulation
----------

First preprocess:

```bash
gmx grompp -f mdp/prd.mdp -po prd -pp prd -o prd
```

The rerun the previous water simulation with the test particle. You must specify
the water simulation's trajectory file with the `-rerun` flag. 

```bash
gmx mdrun -rerun prd.xtc -deffnm prd
```

Analysis
--------


Summary
-------


[Return to main page](https://wesbarnett.github.io/gromacs-tutorials/)
