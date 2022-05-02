# capstone-code

This repo contains matlab code regarding population growth problems. Within the folder `symmetric_pde_model` are four files. 

`basic_sim.m` is a simple model where particle sizes are stored in an array, updated at each time step, and then checked to see if they
are above a threshold size. If so, they split.

`advection_uneven_2D_grid.m` is currently incomplete. It attempts to solve a PDE for a population in which cells require two nutrients
to be above a threshold level before division. It uses an uneven grid in the nutrient levels to deal with splitting accurately.

`simple_advection_1D.m` is the code which generated the figures used in my capstone paper. It successfully solves a simpler PDE
where particles have a single nutrient which must be collected before division. It is, like the rest of the models in this file,
a symmetric method.

`uneven_grid_testing.m` is code for generating and visualizing an uneven grid to assist in completing `advection_uneven_2D_grid.m`.
