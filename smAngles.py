# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 09:24:32 2018
@author: Sula Mueller

DESCRIPTION: calculates all collimator leaf angles for each leaf
detailed description: see MAIN ROUTINE determineAllAnglesAndCoordinates
"""
import numpy as np
import math as m

''' determineX
determine position (coordinate) of (left) edge of current alpha-leaf
i/j: current tile/ pixel
other input: scalar, only value for x-direction '''
def determineX(EAD, edge_Module, edge_Element, delta_e, delta_p, i, j):
    # starts at left edge of first pixel, then considers all left edges
    # right edge of last pixel considered since loop til N+1
    return EAD + edge_Module + edge_Element + i*delta_e + j*delta_p

''' determineZ
determine position (coordinate) of (front) edge of current rho-leaf
k/m/n: current tile/ pixel/ module
other input: scalar, only value for x-direction
difference to determineX: can also consider different rows of modules '''
def determineZ(EAD, size_M, gap_M, edge_Module, edge_Element, delta_e,
               delta_p, k, m, n):  # n... index of module
    return (EAD + n*(size_M+gap_M) + edge_Module + edge_Element + k*delta_e
            + m*delta_p)

''' determineAngle
calculate angle of leaf
INPUT: given distance of current edge to centerline (center of module)
OUTPUT: angle (can be alpha or rho) '''
def determineAngle(b, SDD):
    if b==0: return m.pi/2  # for 90
    return np.arctan(SDD/b)

''' determineAnglesSparseInfo
does the same as determineAllAnglesAndCoordinates with simplifications:
    * edges, gaps = 0
    * isotropic: N_E,P [x] = N_E,P [z] -> scalar input
    * does NOT consider thickness of grid for coordinate calculations ->
      will be off by th_grid/2 (probably low influence) '''
def determineAnglesSparseInfo(SDD, size_M, N_E, N_P, n):
    return determineAllAnglesAndCoordinates(SDD, [size_M/2, size_M/2],
                                            [size_M, size_M], [0, 0], [0, 0],
                                            [0, 0], [0, 0], [0, 0], [0, 0],
                                            [N_E, N_E], [N_P, N_P], n)

''' MAIN ROUTINE: determineAllAnglesAndCoordinates
# DESCRIPTION:
    calculates angels at specified coordinates:
        edges of pixels, tiles and module
    saves angle and position for each leaf (left, front positions)
    considers leaves of grid -> makes 1 entry more than total number of pixels
                                in row
# INPUTS:
    vectors with arg[0] = arg[x], arg[1] = arg[z]
    M... Module, E... Element, P... Pixel
    SDD... Source detector difference
    EAD... Edge Axis Difference -> the distance from central axis (beam) to
           actual edge of module (given as NEGATIVE value)
    "edge": edge of element <-> edge of next functional unit,
    "gap": distance between elements of same kind
    edge_M/E: _0-> CS origin side, _1 -> other side, positive x,z
        (if not specified, assume _0 -> on CO side)
    for several (n) modules in z: assume, they are appended at positive z
# OUTPUT:
    [0]: angles for leaves in z-direction
    [1]: angles for leaves in x-direction
    [2]: x-coordinates of (front) leaf edges in z-direction
    [3]: z-coordinates of (left) leaf edges in x-direction
         REMARK: coordinates of back and right side leaf edges derived from 
                 known leaf thickness
    [4]: number of pixels on tile (needed later for filling of gaps between
         tiles) '''
def determineAllAnglesAndCoordinates(SDD, EAD, size_M, edge_M0, edge_M1,
                                     gap_M, edge_E0, edge_E1, gap_E,
                                     N_E, N_P, n):
    alphas = np.zeros(int(N_E[0])*(int(N_P[0])+1))  # angles of leaves along z
    X = np.zeros(int(N_E[0])*(int(N_P[0])+1))  # (front) edges of these leaves
    rhos = np.zeros(N_E[1]*(N_P[1]+1))  # angles of leaves in x-direction
    Z = np.zeros(N_E[1]*(N_P[1]+1))  # (left) edges of these leaves
    delta_e = []  # distance between 2 Element center points
    delta_e.append((size_M[0]-edge_M0[0]-edge_M1[0])/N_E[0])  # in x
    delta_e.append((size_M[1]-edge_M0[1]-edge_M1[1])/N_E[1])  # in z
    delta_p = []  # distance between 2 Pixel center points
    delta_p.append((delta_e[0]-gap_E[0]-edge_E0[0]-edge_E1[0])/N_P[0]) # in x
    delta_p.append((delta_e[1]-gap_E[1]-edge_E0[1]-edge_E1[1])/N_P[1]) # in z

    #in x:
    for i in range(0, N_E[0]):  # for each tile
        for j in range(0, N_P[0]+1):  # for each pixel
            x = determineX(EAD[0], edge_M0[0], edge_E0[0], delta_e[0],
                           delta_p[0], i, j)
            # x: position/coordinate of (front) edge of current leaf
            alphas[i*(N_P[0]+1)+j] = determineAngle(x,SDD)
            X[i*(N_P[0]+1)+j] = x
    #in z:
    for i in range(0, N_E[1]):  # for each tile
        for j in range(0,N_P[1]+1):  # for each pixel
            z = determineZ(EAD[1], size_M[1], gap_M, edge_M0[1], edge_E0[1],
                           delta_e[1], delta_p[1], i, j, n)
            # z: position/coordinate of (left) edge of current leaf
            rhos[i*(N_P[1]+1)+j] = determineAngle(z,SDD)
            Z[i*(N_P[1]+1)+j] = z
    return [alphas, rhos, X, Z, N_P]