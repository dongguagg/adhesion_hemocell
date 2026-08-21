# HemoCell RBC Adhesion Extension

This project extends HemoCell with adhesion interactions between red blood cells (RBCs) and between RBCs and solid walls. Both interactions combine short-range Lennard-Jones repulsion with medium-range Morse attraction, while using independent parameters for cell-cell and cell-wall adhesion.

The repository also contains three small example cases:

- `twoCellShear`: two adhered RBCs under Couette shear, with the lower RBC anchored through wall adhesion.
- `twoCellPull`: a pull-off test in which an external force is applied to one of two adhered RBCs.
- `bifurcation`: adhesive RBC flow through a bifurcating vessel, using RBC-RBC adhesion and solid-wall repulsion.

The modified HemoCell source is located in `src_hemocell`. Each example directory contains its own configuration, RBC files, build script, and run script.
