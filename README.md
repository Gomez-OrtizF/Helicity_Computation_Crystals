# HelicityComputation
This script is intended for the computation of the Helicity of crystalline structures.
A full description of the method, their strengths and weaknesses and the application to real systems 
can be found in XXXXXXXXXXXXX(TO BE PUBLISHED)XXXXXXXXXXXXXXXXX 

The script takes three parameters:
    a (int):   The number of atoms in the home unit cell of the materials.
    b (string):The path to a file containing the atomic positions in fractional units
               of the reference high symmetry structure in the space group of the low symmetry phase.
    c (string):The path file containing the atomic positions in fractional units
             of the low symmetry phase in the same order as the previous file (one to one mapping).
And returns the value of the helicity for the crystalline phase
    float: The helicity of the structure.
