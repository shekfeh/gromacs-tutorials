[Return to main page](https://wesbarnett.github.io/gromacs-tutorials/)

Tutorial 6: Methane excess chemical potential through test particle insertion
======================================================

In this tutorial we'll be using test particle insertion (TPI) to calculate the excess
chemical potential of methane solvation in water. Most users are unaware that
GROMACS has a built-in method for running TPI. This tutorial will not be a
comprehensive discussion on the statistical mechanics of TPI, but will address
issues when needed. The user is encouraged to seek out scientific resources
regarding this method

TPI involves perturbing some state to some other, very similar state. We will be
taking bulk water and attempting to insert a methane particle and measuring the
potential energy change from this. So for us, state A is the bulk water system,
and state B is the water system with a methane. 

With GROMACS you need to run state A as a normal MD simulation. We already did
this for our case of bulk water in Tutorial 1. We'll reuse the output trajectory
files for inserting the methane.

Setup
-----

### Create water system

Follow Tutorial 1 to run a system containing TIP4PEW water.

### Add test particle to topology file

Our original topology file just had water. In the new topology file we simply
need to add 1 test particle, and it needs to be the last molecule in the system.
We'll use `opls_066` for the particle's atom type which is OPLS's united atom
methane. Here's what my final topology file looks like (the number of waters
will be different for your system):

    #include "oplsaa.ff/forcefield.itp"
    #include "oplsaa.ff/tip4pew.itp"

    [ moleculetype ]
    ; Name          nrexcl
    Methane         3

    [ atoms ]
    ;   nr       type      resnr residue  atom   cgnr     charge       mass
         1       opls_066  1     CH4      C      1          0     16.043


    [ System ]
    TIP4PEW in water

    [ Molecules ] 
    SOL               395
    Methane           1

### Add test particle to gro file

You also need to add the test particle to the gro file. Simply edit `conf.gro`
(or any of the other `.gro` files uses) and add a line at the end containing the
test particle's position (right before the box coordinates). The line I added
looks like this:

    396CH4      C 1581   0.000   0.000   0.000

Additionally you need to add 1 to the total number of particles in the system on
the second line of the `.gro` file.

### Parameter files

We only need one parameter file for TPI. Simply copy `prd.mdp` from your bulk
water simulation and change `integrator` to `tpi`. You should change `nsteps` to
the number of insertions per frame that you want to attempt. You will also need
to change `cutoff-scheme` to `group`, since `Verlet` has not be implemented for
TPI.

Simulation
----------

For the simulation we are just rerunning the bulk water simulation using the
saved trajectory file (which was named `prd.xtc` in the first tutorial). To do
this first run `grompp`:

```bash
$ gmx grompp -f mdp/tpi.mdp -o tpi -po tpi -pp tpi -c conf.gro
```

Now use the `-rerun` flag with `mdrun`:

```bash
$ gmx mdrun -deffnm tpi -rerun prd.xtc
```

Analysis
--------

The log file file should contain a l

Summary
-------


[Return to main page](https://wesbarnett.github.io/gromacs-tutorials/)
