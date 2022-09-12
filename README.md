# Simulation of the model introduced in Macroscopic and microscopic dynamics of a pedestrian cross flow: Part-II modelling

WORK IN PROGRESS!

simulator\_cir is the simulator for circular models, simulator_ell for elliptical models. Currently, it's basically the same code (disks are defined as ellipses with A=B), but we plan to upload an optimised code for circular models.

simulatorr\_cir uses the solution best\_cir, simulator\_ell uses best\_ell. Different solutions can be copied from the solutions folder.

graphictools includes graphics functions, group pedestrian behaviour functions (it includes group behaviour functions, but they are not used) and ellipses includes physical dynamics functions.

SDL2 is used for graphics

To compile

g++ -o simulator_* simulator_*.cpp graphictools.cpp `sdl2-config --cflags --libs`




