# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 08:31:13 2018
@author: Sula Mueller

DESCRIPTION: turns information about leaf angles into Facets for STLA File
# each facet is defined by: normal, vertex1, vertex2, vertex3
    * normal: nx, ny, nz; always defined in POSITIVE direction
      (routine flips them automatically for faces facing opposite)
    * vertex: x,y,z
    => 12 values for 1 triangle
# need 2 triangles to define rectangle
# all values are stored in linear fashion -> [nx,ny,nz,v0,v1...,nx,...]
# leaves orientated along x-axis ("rho-leaves") are continuous, "alpha-leaves"
    are segmented to avoid overleafing
"""
import math as m

''' getD
* leaves are tilted -> upper edges are shifted by dx, dz
* getD to calculate that shift
* works in x,z alike '''
def getD(angle, h_grid):
    if angle == m.pi/2:
        dx = 0
    else:
        dx = h_grid/ m.tan(angle)
    return -dx

''' appends values in vector v to existing array F '''
def appendVector(F, v):
    s = len(v)
    for i in range(0, s):
        F.append(v[i])
    return F

''' flips 3-dimensional vector to turn into its negative '''
def neg(n):
    return [-n[0], -n[1], -n[2]]

''' appends F by values for two triangles forming a rectangle '''
def appendRectangle(F, n, x0, x1, x2, x3):
    F = appendVector(F, n)
    F = appendVector(F, x0)
    F = appendVector(F, x1)
    F = appendVector(F, x2)
    
    F = appendVector(F, n)
    F = appendVector(F, x1)
    F = appendVector(F, x3)
    F = appendVector(F, x2)
    return F

''' returns 6x2 surfaces of "cube" defined by 8 points and 3 normals '''
def CubeFromPoints(n_a, n_y, n_b, x0, x1, x2, x3, x0_, x1_, x2_, x3_):
    F = []
    F = appendRectangle(F, neg(n_a), x0, x1, x2, x3)
    F = appendRectangle(F, n_a, x0_, x1_, x2_, x3_)
    F = appendRectangle(F, neg(n_y), x0, x0_, x1, x1_)
    F = appendRectangle(F, n_y, x2, x2_, x3, x3_)
    F = appendRectangle(F, neg(n_b), x0, x0_, x2, x2_)
    F = appendRectangle(F, n_b, x1, x1_, x3, x3_)
    return F

''' same as CubeFromPoints with one pair of faces not parallel (b)
    -> need one more normal '''
def PolyederFromPoints(n_a, n_y, n_b0, n_b1, x0, x1, x2, x3, x0_, x1_, x2_,
                       x3_):
    return irregPolyederFromPoints(n_a, n_a, n_y, n_b0, n_b1, x0, x1, x2, x3,
                                   x0_, x1_, x2_, x3_)

''' irregPolyederFromPoints
* returns 6x2 surfaces of irregular "cube" defined by 8 points and 5 normals
* 2 pairs of faces NOT parallel, third pair IS parallel (top, bottom)
* a,b: if alpha-> a=x, b=z
       if rho  -> a=z, b=x '''
def irregPolyederFromPoints(n_a0, n_a1, n_y, n_b0, n_b1, x0, x1, x2, x3, x0_,
                            x1_, x2_, x3_ ):
    F = []
    F = appendRectangle(F, neg(n_a0), x0, x1, x2, x3)
    F = appendRectangle(F, n_a1, x0_, x1_, x2_, x3_)
    F = appendRectangle(F, neg(n_y), x0, x0_, x1, x1_)
    F = appendRectangle(F, n_y, x2, x2_, x3, x3_)
    F = appendRectangle(F, neg(n_b0), x0, x0_, x2, x2_)
    F = appendRectangle(F, n_b1, x1, x1_, x3, x3_)
    return F

''' FacetizeOneLeafRho
DESCRIPTION:
    returns facet values of leaf orientated along x-axis (rho-orientation)
INPUT:
    dx0 -> tilt of left (small) edge face
    dx1 -> tilt of right (small) edge face
    dz -> tilt of entire leaf (difference between bottom, top coordinates)
    rho -> angle of tilt (causing dz)
    X_0 -> [X[0],Z[i]] -> left front bottom point (i... leaf)
    X_1 -> [X[end]+th_grid,Z[i]] -> right front bottom point '''
def FacetizeOneLeafRho(dx0, dx1, dz, rho, X_0, X_1, h_grid, th_grid):
    ''' determine normals: '''
    abs_n0 = m.sqrt(h_grid*h_grid + dx0*dx0)  # length of normals...
    abs_n1 = m.sqrt(h_grid*h_grid + dx1*dx1)  # ...for normalization
    # normalized normal of left (small) edge face:
    n_x0 = [h_grid/abs_n0, -dx0/abs_n0, 0]
    # normalized normal of right (small) edge face:
    n_x1 = [h_grid/abs_n1, -dx1/abs_n1, 0]  # pointing to positive z
    n_y = [0, 1, 0]
    n_z = [0, -m.cos(rho), m.sin(rho)]  # pointing to positive z
    ''' determine 8 edge points of leaf: '''
    x0 = [X_0[0], 0, X_0[1]]
    x1 = [X_1[0], 0, X_1[1]]
    x2 = [X_0[0]+dx0, h_grid, X_0[1]+dz]
    x3 = [X_1[0]+dx1, h_grid, X_1[1]+dz]
    # +th_grid: cuts them short in z to avoid overleafing with rho-leaves
    x0_ = [X_0[0], 0, X_0[1]+th_grid]
    x1_ = [X_1[0], 0, X_1[1]+th_grid]
    x2_ = [X_0[0]+dx0, h_grid, X_0[1]+dz+th_grid]
    x3_ = [X_1[0]+dx1, h_grid, X_1[1]+dz+th_grid]
    return PolyederFromPoints(n_z, n_y, n_x0, n_x1, x0, x1, x2, x3, x0_, x1_,
                              x2_, x3_)

''' FacetizeOneLeafAlpha
DESCRIPTION:
    returns facet values of leaf oriented along z-axis (alpha-orientation)
INPUT:
    dz0 -> tilt of front (small) edge face
    dz1 -> tilt of back (small) edge face
    alpha -> angle of tilt (causing dx)
    X_0 -> [X[j], Z[i]+th_grid] -> left front bottom point
           Z[i]+th_grid: cuts them short in z to avoid overleafing with
                         rho-leaves
    X_1 -> [X[j], Z[i+1]] -> left back bottom point
REMARK: i -> row (corresponding rho leaf), j -> column '''
def FacetizeOneLeafAlpha(dz0, dz1, alpha, X_0, X_1, h_grid, th_grid):
    # tilt of entire leaf (difference between bottom, top coordinates):
    dx = getD(alpha, h_grid)
    ''' determine normals: '''
    n_x = [m.sin(alpha), -m.cos(alpha), 0]  # pointing to positive x
    n_y = [0, 1, 0]
    abs_n0 = m.sqrt(h_grid*h_grid + dz0*dz0)  # length of normals...
    abs_n1 = m.sqrt(h_grid*h_grid + dz1*dz1)  # ...for normalization
    # normalized normal of front (small) edge face:
    n_z0 = [0, -dz0/abs_n0, h_grid/abs_n0]
    # normalized normal of back (small) edge face:
    n_z1 = [0, -dz1/abs_n1, h_grid/abs_n1]
    ''' determine 8 edge points of leaf: '''
    x0 = [X_0[0], 0, X_0[1]]
    x1 = [X_1[0], 0, X_1[1]]
    x2 = [X_0[0]+dx, h_grid, X_0[1]+dz0]
    x3 = [X_1[0]+dx, h_grid, X_1[1]+dz1]
    # +th_grid: cuts them short in z to avoid overleafing with rho-leaves
    x0_ = [X_0[0]+th_grid, 0, X_0[1]]
    x1_ = [X_1[0]+th_grid, 0, X_1[1]]
    x2_ = [X_0[0]+th_grid+dx, h_grid, X_0[1]+dz0]
    x3_ = [X_1[0]+th_grid+dx, h_grid, X_1[1]+dz1]
    return PolyederFromPoints(n_x, n_y, n_z0, n_z1, x0, x1, x2, x3, x0_, x1_,
                              x2_, x3_)

''' FacetizeGapRho
DESCRIPTION:
    to totally fill the gap between 2 tiles with one (thick) leaf
    melts last leaf of previous tile and first leaf of next tile
    returns Facet vector for that structure
INPUT:
    dx0 -> tilt of left (small) edge face
    dx1 -> tilt of right (small) edge face
    dz -> tilt of front face ("current leaf")
    dz1 -> tilt of back face ("next leaf" (-> i+1))
    rho/ rho1 -> angles of tilt (causing dzs) of front/back faces
    X_0 -> [X[0],Z[i]] -> left front bottom point (i... "current leaf")
    X_1 -> [X[end]+th_grid, Z[i]] -> right front bottom point
    X_4 -> [X[0], Z[i+1]+th_grid] -> left back bottom point ("next leaf") '''
def FacetizeGapRho(dx0, dx1, dz, dz1, rho, rho1, X_0, X_1, X_4_, h_grid):
    ''' determine normals: '''
    abs_n0 = m.sqrt(h_grid*h_grid + dx0*dx0)  # length of normals...
    abs_n1 = m.sqrt(h_grid*h_grid + dx1*dx1)  # ...for normalization
    # normalized normal of left (small) edge face:
    n_x0 = [h_grid/abs_n0, -dx0/abs_n0, 0]
    # normalized normal of right (small) edge face:
    n_x1 = [h_grid/abs_n1, -dx1/abs_n1, 0]
    n_y = [0, 1, 0]
    n_z0 = [0, -m.cos(rho), m.sin(rho)]  # pointing to positive z
    n_z1 = [0, -m.cos(rho1), m.sin(rho1)]  # pointing to positive z
    ''' determine 8 edge points of leaf: '''
    x0 = [X_0[0], 0, X_0[1]]
    x1 = [X_1[0], 0, X_1[1]]
    x2 = [X_0[0]+dx0, h_grid, X_0[1]+dz]
    x3 = [X_1[0]+dx1, h_grid, X_1[1]+dz]
    # +th_grid: cuts them short in z to avoid overleafing with rho-leaves
    x4_ = [X_0[0], 0, X_4_[1]]
    x5_ = [X_1[0], 0, X_4_[1]]
    x6_ = [X_0[0]+dx0, h_grid, X_4_[1]+dz1]
    x7_ = [X_1[0]+dx1, h_grid, X_4_[1]+dz1]
    return irregPolyederFromPoints(n_z0, n_z1 ,n_y, n_x0, n_x1, x0, x1, x2, x3,
                                   x4_, x5_, x6_, x7_)

''' FacetizeGapAlpha
DESCRIPTION:
    to totally fill the gap between 2 tiles with one (thick) leaf
    melts last leaf of previous tile and first leaf of next tile
    returns Facet vector for that structure
INPUT:
    dz0 -> tilt of front (small) edge face
    dz1 -> tilt of back (small) edge face
    alpha/ alpha1 -> angles of tilt (causing dx) of right/ left faces
    X_0 -> [X[j], Z[i]+th_grid] -> left front bottom point (j -> "current 
           leaf")
    X_1 -> [X[j], Z[i+1]] -> left back bottom point
    X_4 -> [X[j+1]+th_grid, Z[i]+th_grid] -> right front bottom point (of "next
           leaf")
    i... row (corresponding rho leaf), j... column '''
def FacetizeGapAlpha(dz0, dz1, alpha, alpha1, X_0, X_1, X_4_, h_grid):
    # tilt of left face (difference between bottom, top coordinates):
    dx = getD(alpha, h_grid)
    # tilt of right face (difference between bottom, top coordinates):
    dx1 = getD(alpha1, h_grid)
    ''' determine normals: '''
    n_x0 = [m.sin(alpha), -m.cos(alpha), 0]  # pointing to positive x
    n_x1 = [m.sin(alpha1), -m.cos(alpha1), 0]  # pointing to positive x
    n_y = [0, 1 ,0]
    abs_n0 = m.sqrt(h_grid*h_grid + dz0*dz0)  # length of normals...
    abs_n1 = m.sqrt(h_grid*h_grid + dz1*dz1)  # ...for normalization
    # normalized normal of front (small) edge face:
    n_z0 = [0, -dz0/abs_n0, h_grid/abs_n0]
    # normalized normal of back (small) edge face:
    n_z1 = [0, -dz1/abs_n1, h_grid/abs_n1]
    ''' determine 8 edge points of leaf: '''
    x0 = [X_0[0], 0, X_0[1]]
    x1 = [X_1[0], 0, X_1[1]]
    x2 = [X_0[0]+dx, h_grid, X_0[1]+dz0]
    x3 = [X_1[0]+dx, h_grid, X_1[1]+dz1]
    # +th_grid: cuts them short in z to avoid overleafing with rho-leaves
    x4_ = [X_4_[0], 0, X_0[1]]
    x5_ = [X_4_[0], 0, X_1[1]]
    x6_ = [X_4_[0]+dx1, h_grid, X_0[1]+dz0]
    x7_ = [X_4_[0]+dx1, h_grid, X_1[1]+dz1]
    return irregPolyederFromPoints(n_x0, n_x1, n_y, n_z0, n_z1, x0, x1, x2, x3,
                                   x4_, x5_, x6_, x7_)

''' FacetizeEntireGrid
DESCRIPTION:
    turns information about leaf angles into Facets for STLA File
    returns entire Facets vector of ALL facets in linear fashion
    -> [nx, ny ,nz, v0, v1..., nx,...]
INPUT:
    Angles (from smAngles.determineAllAnglesAndCoordinates)
    Grid Ratio, thickness of grid leaves '''
def FacetizeEntireGrid(Angles, GR, th_grid):
    ''' INITIALIZE: '''
    alphas = Angles[0]  # per output of determineAllAnglesAndCoordinates
    rhos = Angles[1]
    X = Angles[2]
    Z = Angles[3]
    N_P = Angles[4]
    
    a = len(alphas)
    r = len(rhos)
    
    h_grid = th_grid*GR
    F = []  # Facets vector with all values
    
    ''' EDGE RECYCLING: reuse tilt of leaves '''
    # tilt of first/last alpha leaf determines tilt of small edge faces
    # of rho-leaves (dx0, dx1):
    dx0 = getD(alphas[0], h_grid)
    dx1 = getD(alphas[a-1], h_grid)
    # tilt of "current" rho-leaf determines tilt of small front edge face
    # of alpha-leaves (dz0):
    dz0 = getD(rhos[0], h_grid)
    # REMARK: dz1 will be calculated FIRST for back edge face of alpha-leaves,
    # then reused for next rho-leaf with dz0 = dz1 at end of alpha-loop
    skip_Flag_rho = False
    # gapfiller already includes next leaf, next leaf will be skipped
    
    ''' RHO LEAVES: '''
    for i in range(0, r):  # z
        if skip_Flag_rho:  # skip all leaves that are included in "gap-filler"
            skip_Flag_rho = False  # continue with leaf after that
        else:
            rho = rhos[i]  # from input Angles
            if ((i+1)%(N_P[1]+1)==0 and (i<r-1) ):  # gap, but not last leaf
                rho1 = rhos[i+1]  # of "next" leaf (other side of gap)
                dz1 = getD(rho1, h_grid)
                F = appendVector(F, FacetizeGapRho(dx0, dx1, dz0, dz1, rho,
                                                   rho1, [X[0], Z[i]],
                                                   [X[a-1]+th_grid, Z[i]],
                                                   [X[0], Z[i+1]+th_grid],
                                                   h_grid))
                dz0 = dz1  # for use in alpha leaves and reuse next rho leaf
                skip_Flag_rho = True  # to skip next leaf
                # (because already included in edge)
            else:  # normal leaf
                F = appendVector(F, FacetizeOneLeafRho(dx0, dx1, dz0, rho,
                                                       [X[0], Z[i]],
                                                       [X[a-1]+th_grid, Z[i]],
                                                       h_grid, th_grid))
                if i==r-1: return F  # return here, because last rho-leaf
                                     # can't be boardered by more alpha-leaves
            ''' ALPHA LEAVES: included in rho loop because dont need them if
                gap filled '''
            if skip_Flag_rho:
                i = i+1 # so that alpha leave boarder to second "leaf" of gap
            dz1 = getD(rhos[i+1], h_grid)  # needed for back (small) edge face
                                           # reuse it later for rho
            skip_Flag_alpha = False  # x-direction has tile gaps too
            for j in range(0, a):  # x      
                alpha = alphas[j]
                if (  (j+1)%(N_P[0]+1)==0 and (j<a-1)):  # gap
                    alpha1 = alphas[j+1]
                    # of "next" leaf; other side of gap
                    F = appendVector(F, FacetizeGapAlpha(dz0, dz1, alpha,
                                                         alpha1,
                                                         [X[j], Z[i]+th_grid],
                                                         [X[j], Z[i+1]],
                                                         [X[j+1]+th_grid,
                                                          Z[i]+th_grid],
                                                         h_grid))
                    skip_Flag_alpha = True
                else:
                    if skip_Flag_alpha:  # skip next leaf
                        skip_Flag_alpha = False
                    else:
                        F = appendVector(F, FacetizeOneLeafAlpha(dz0, dz1,
                                                                 alpha,
                                                                 [X[j], 
                                                                 Z[i]+th_grid],
                                                                 [X[j],Z[i+1]],
                                                                 h_grid,
                                                                 th_grid))
            dz0 = dz1  # for reuse in rho leaves