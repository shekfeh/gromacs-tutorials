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

### Parameter files

We only need on parameter file for TPI.
TODO: expand

TODO: update
Here's an explanation of the new parameters that are used in each file:

|  parameter     | value     | explanation |
| ---------------|-----------|-------------| 
|  pull          | yes | Use the pull code. |
|  pull-ngroups       | 2         | We have two groups that we're pulling.  |
|  pull-group1-name   | CA     | We specified this in the index file. For us this will be the carbon of one of the methanes, although we probably could have chosen the entire methane. If we did that it would have been pulled along the COM of the entire molecule. |
|  pull-group2-name   | CB  | The carbon of the other methane. |
|  pull-ncoords       |  1  | We are pulling along only one coordinate. |
|  pull-coord1-geometry  | distance  | We're going to pull along the vector connecting our two groups. |
| pull-coord1-type    |   umbrella | Use an umbrella (harmonic) potential for this coordinate.| 
|  pull-coord1-groups | 1 2       | For this pull coordinate these are the two groups (defined below) which will be pulled. You can actually have more thane one pull coordinate and so do pulling across different sets of molecules, but that's not applicable here. |
|  pull-coord1-k      | 5000.0  | The force constant used in the umbrella potential in kJ/(mol nm). |
|  pull-coord1-init   | WINDOW  | This is the distance we want our two groups to be apart. I've put this keyword here that I'll replace in our bash script for each window |
|  pull-coord1-rate   | 0.0    | We don't want the groups to move along the coordinate any, so this is 0. |
|  pull-coord1-start         | no        | We're manually specifying the distance for each window, so we do not want to add the center of mass distance to the calculation. |

### Topology file



Simulation
----------

TODO: make sure these are correct
For the simulation we are just rerunning the bulk water simulation using the
saved trajectory file. To do this first run `grompp`:

```bash
$ gmx grompp -f tpi.mdp -o tpi -po tpi -pp tpi -c tpi.gro
```

Now use the `-rerun` flag with `mdrun`:

```bash
$ gmx mdrun -deffnm tpi -rerun
```

Analysis
--------

The log file file should contain a line with Greek letter mu.
TODO: expand

Summary
-------


[Return to main page](https://wesbarnett.github.io/gromacs-tutorials/)
