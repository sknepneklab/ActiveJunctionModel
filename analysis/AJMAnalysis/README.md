# This is a toolset for analysing results of AJM simulations

This module contains several submodules

# Grander module

This summodule is a part of the Active Junction Model analysis toolkit. Its purpose is to 
implement analysis of Vertex Model simulations closely following this paper:

F. Graner, B. Dollet, C. Raufaste and P. Marmottant, Discrete rearranging disordered patterns, part I:
Robust statistical tools in two or three dimensions, Eur. Phys. J. E 25, 349–369 (2008).

# utils/HalfEdge

This is a utility module that read in mesh from a JSON file and stores it as half-edge data structure. 
It is used by several submodules in this analysis toolkit.