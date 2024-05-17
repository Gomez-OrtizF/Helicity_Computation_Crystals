import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import math
def comp_hel(a, b, c):
    """
    This function takes three parameters.

    Parameters:
    a (int):   The number of atoms in the home unit cell of the materials.
    b (string):The path to a file containing the atomic positions in fractional units
               of the reference high symmetry structure in the space group of the low symmetry phase.
    c (string):The path file containing the atomic positions in fractional units
             of the low symmetry phase in the same order as the previous file (one to one mapping).


    And returns the value of the helicity 
    Returns:
    float: The helicity of the structure.
    """
    #-------------------------------------------------------------------
    #               INITIALIZATION
    #-------------------------------------------------------------------
    N_at=a
    coords_ch1=np.zeros((N_at,3))
    ref_hs=np.zeros((N_at,3))
    vel_1=np.zeros((N_at,3))
    helicity_1=0.0
    with open(b,'r') as f:
        for i in range(N_at):
            line=f.readline().split()
            ref_hs[i][0]=float(line[0])
            ref_hs[i][1]=float(line[1])
            ref_hs[i][2]=float(line[2])
    with open(c,'r') as f:
        for i in range(N_at):
            line=f.readline().split()
            coords_ch1[i][0]=float(line[0])
            coords_ch1[i][1]=float(line[1])
            coords_ch1[i][2]=float(line[2])
            vel_1[i][0]=coords_ch1[i][0]-ref_hs[i][0]
            vel_1[i][1]=coords_ch1[i][1]-ref_hs[i][1]
            vel_1[i][2]=coords_ch1[i][2]-ref_hs[i][2]
            for j in enumerate(vel_1[i]):
                if j[1]>0.5:
                    vel_1[i][j[0]]=vel_1[i][j[0]]-1
                elif j[1]<-0.5:
                    vel_1[i][j[0]]=vel_1[i][j[0]]+1
    neig=np.zeros((N_at,6))
    vec1=np.zeros(3)
    vec2=np.zeros(3) 
    vec3=np.zeros(3)
    #-------------------------------------------------------------------
    #               COMPUTATION
    #-------------------------------------------------------------------
    neig=np.zeros((N_at,6))
    vec1=np.zeros(3)
    vec2=np.zeros(3)
    vec3=np.zeros(3)
    for i in range(N_at):
        #We start looking for neighbors
        diff=np.zeros(3)
        derivative1=np.zeros(3)
        derivative2=np.zeros(3)
        derivative3=np.zeros(3)
        dist1=1.0
        for j in range(N_at):
            diff[0]=ref_hs[j][0]-ref_hs[i][0]
            diff[1]=ref_hs[j][1]-ref_hs[i][1]
            diff[2]=ref_hs[j][2]-ref_hs[i][2]
            for k in enumerate(diff):
                if k[1]>0.5:
                    diff[k[0]]=diff[k[0]]-1
                elif k[1]<-0.5:
                    diff[k[0]]=diff[k[0]]+1
            if (np.linalg.norm(diff)<dist1 and i!=j):
                dist1=np.linalg.norm(diff)
                vec1[:]=diff[:]
                neig[i][0]=j
                neig[i][1]=dist1
        #We already have the first neighbors of each atom now lets find now
        #the 2nd nearest neighbor in a different direction
        diff2=np.zeros(3)
        dist2=1.0
        for j in range(N_at):
            diff2[0]=ref_hs[j][0]-ref_hs[i][0]
            diff2[1]=ref_hs[j][1]-ref_hs[i][1]
            diff2[2]=ref_hs[j][2]-ref_hs[i][2]
            for k in enumerate(diff2):
                if k[1]>0.5:
                    diff2[k[0]]=diff2[k[0]]-1
                elif k[1]<-0.5:
                    diff2[k[0]]=diff2[k[0]]+1
            matrix=np.zeros((2,3))
            matrix[0][:]=vec1[:]
            matrix[1][:]=diff2[:]
            if (np.linalg.norm(diff2)<dist2 and i!=j and np.linalg.matrix_rank(matrix)==2):
                dist2=np.linalg.norm(diff2)
                neig[i][2]=j
                neig[i][3]=dist2
                vec2[:]=diff2[:]
        #We already have the second neighbors of each atom now lets find
        #the 3rd nearest neighbor in a different direction
        diff3=np.zeros(3)
        dist3=1.0
        for j in range(N_at):
            diff3[0]=ref_hs[j][0]-ref_hs[i][0]
            diff3[1]=ref_hs[j][1]-ref_hs[i][1]
            diff3[2]=ref_hs[j][2]-ref_hs[i][2]
            for k in enumerate(diff3):
                if k[1]>0.5:
                    diff3[k[0]]=diff3[k[0]]-1
                elif k[1]<-0.5:
                    diff3[k[0]]=diff3[k[0]]+1
            matrix=np.zeros((3,3))
            matrix[0][:]=vec1[:]
            matrix[1][:]=vec2[:]
            matrix[2][:]=diff3[:]
            if (np.linalg.norm(diff3)<dist3 and i!=j and np.linalg.matrix_rank(matrix)==3):
                dist3=np.linalg.norm(diff3)
                neig[i][4]=j
                neig[i][5]=dist3
                vec3[:]=diff3[:]
        #-------------------------------#
        #Now we have for each atom a set of 3 independent directions
        #we can compute derivatives along those directions and extract    
        #the values of the derivatives along cartesian directions: x,y,z
        #-------------------------------#
        #we dont normalize because the length of the vector is important for solving linear system
        derivative1[:]=(vel_1[int(neig[i][0])][:]-vel_1[i][:])
        derivative2[:]=(vel_1[int(neig[i][2])][:]-vel_1[i][:])
        derivative3[:]=(vel_1[int(neig[i][4])][:]-vel_1[i][:])
        #-------------------------------#
        #Now we have the values of the directional derivatives
        #In order to compute the values along cartesian directions
        #we need to solve the linear system of equations
        #-------------------------------#
        parx_x=0.0
        parx_y=0.0
        parx_z=0.0
        pary_x=0.0
        pary_y=0.0
        pary_z=0.0
        parz_x=0.0
        parz_y=0.0
        parz_z=0.0
        coef=np.zeros((3,3))
        coef[0][:]=vec1[:]
        coef[1][:]=vec2[:]
        coef[2][:]=vec3[:]
        temp=np.linalg.solve(coef,np.array([derivative1[0],derivative2[0],derivative3[0]]))
        parx_x=temp[0]
        parx_y=temp[1]
        parx_z=temp[2]
        temp=np.linalg.solve(coef,np.array([derivative1[1],derivative2[1],derivative3[1]]))
        pary_x=temp[0]
        pary_y=temp[1]
        pary_z=temp[2]
        temp=np.linalg.solve(coef,np.array([derivative1[2],derivative2[2],derivative3[2]]))
        parz_x=temp[0]
        parz_y=temp[1]
        parz_z=temp[2]
        #With the values of the partial derivatives we compute the helicity
        helicity_1+=vel_1[i][0]*(parz_y-pary_z)+vel_1[i][1]*(parx_z-parz_x)+vel_1[i][2]*(pary_x-parx_y)
    return helicity_1/N_at
def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Calculate the helicity of a given structure.")
    # Add arguments
    parser.add_argument('a', type=int, help="The number of atoms")
    parser.add_argument('b', type=str, help="The path to fractional coordinates of reference in the low symmetry phase")
    parser.add_argument('c', type=str, help="The path to fractional coordinates of structure")
    # Parse the arguments
    args = parser.parse_args()
    # Calculate the sum
    result = comp_hel(args.a, args.b, args.c)
    # Print the result
    print(f"The Helicity of the structure is {result}")
if __name__ == "__main__":
    main()

