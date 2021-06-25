# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 12:00:03 2018

DESCRIPTION: turns information about leaf angles into Lines and Boundaries for 
    IGES File
INPUT: Angles from smAngles
OUTPUT: Matrix -> NxM = (N_P_total+1)*12 x 6
    each entry [:,1:6] represents one line [:,[x0,y0,z0,x1,y1,z1]]
    12 lines for each pixel: 4 bottom lines, 4 vertical lines, 4 top lines
    12 lines for outer edges

@author: Sula Mueller
"""


import numpy as np
import math as m 

''' getD
* leaves are tilted -> upper edges are shifted by dx, dz
* getD to calculate that shift 
* works in x,z alike '''
def getD (angle, h_grid):
    if angle == m.pi/2: return 0 # for 90 -> perpendicular; no tilt
    return h_grid/ m.tan(angle) 

'''getALLPoints
DESCRIPTION: returns line array of all points surrounding pixels
INPUT: coordinates of 4 bottom points ([x0,z0], [x0,z1], [x1,z1], [x1,z0]) 
             x0, z0 should be +th_grid
       dx, dz -> tilt of upper points
OUTPUT: line array '''
def getPixelPoints (x0,x1,z0,z1, dx0,dx1,dz0,dz1, h_grid):
    P = np.zeros((8,3))
    P[0,:] = [x0,0,z0]
    P[1,:] = [x1,0,z0]
    P[2,:] = [x1,0,z1]
    P[3,:] = [x0,0,z1]
    P[4,:] = [x0+dx0,h_grid,z0+dz0]
    P[5,:] = [x1+dx1,h_grid,z0+dz0]
    P[6,:] = [x1+dx1,h_grid,z1+dz1]
    P[7,:] = [x0+dx0,h_grid,z1+dz1]
    return P

'''getPixelLines
DESCRIPTION: returns line array of all lines surrounding pixels
INPUT: coordinates of 4 bottom points ([x0,z0], [x0,z1], [x1,z1], [x1,z0]) 
             x0, z0 should be +th_grid
       dx, dz -> tilt of upper points
OUTPUT: line array
REMARK: DON'T do for last alpha,rho (extra routine for outer lines)'''
def getPixelLines (x0,x1,z0,z1, dx0,dx1,dz0,dz1, h_grid):
    L = np.zeros((12,6))
    '''BOTTOM:'''
    L[0,:] = [x0,0,z0,  x1,0,z0]
    L[1,:] = [x1,0,z0,  x1,0,z1]
    L[2,:] = [x1,0,z1,  x0,0,z1]
    L[3,:] = [x0,0,z1,  x0,0,z0]
    '''VERTICAL:'''
    L[4,:] = [x0,0,z0,  x0+dx0,h_grid,z0+dz0]
    L[5,:] = [x1,0,z0,  x1+dx1,h_grid,z0+dz0]
    L[6,:] = [x1,0,z1,  x1+dx1,h_grid,z1+dz1]
    L[7,:] = [x0,0,z1,  x0+dx0,h_grid,z1+dz1]
    '''TOP:'''
    L[8,:] = [x0+dx0,h_grid,z0+dz0,  x1+dx1,h_grid,z0+dz0]
    L[9,:] = [x1+dx1,h_grid,z0+dz0,  x1+dx1,h_grid,z1+dz1]
    L[10,:] = [x1+dx1,h_grid,z1+dz1,  x0+dx0,h_grid,z1+dz1]
    L[11,:] = [x0+dx0,h_grid,z1+dz1,  x0+dx0,h_grid,z0+dz0]
    return L    
     
''' getAllLinesLoop
DESCRIPTION: turns information about leaf angles into Lines for IGES File
INPUT: Angles from smAngles
OUTPUT Lines: Matrix [(N_P_total+1)*12, 6] 
    each entry [:, 1:6] represents one line [:,[x0,y0,z0,x1,y1,z1]]
    12 lines for each pixel: 4 bottom lines, 4 vertical lines, 4 top lines
    all lines linearly stored (12 entries belong together)      
inner vs outer:
    inner: lines directly boardering the pixels
    outer: edges of entire grid (stored last in same array) 
def getAllLinesLoop (Angles, GR, th_grid):    
    #INITIALIZE:
    alphas = Angles[0]  # per output of determineAllAnglesAndCoordinates
    rhos = Angles[1]
    X = Angles[2]
    Z = Angles[3]
    N_P = Angles[4]
    
    a = len(alphas)
    r = len(rhos)
    
    h_grid = th_grid*GR
    
    # N_P+1 angles on tile -> (a-1)/ (N_P[0] +1) = N_Tile
    N_P_x = ((a-1)/ (N_P[0] +1)) * N_P[0] 
    N_P_z = ((r-1)/ (N_P[1] +1)) * N_P[1]
    N_P_total = N_P_x * N_P_z
    
    L = np.zeros[(N_P_total+1)*12,6] # Lines array
    index = 0 # to count entries in L (despite gaps)
                 
    #INNER LINES:
    dx00 = getD(alphas[0], h_grid) # tilt of left outer edge
    dx0 = dx00                     # tilt of first leaf = tilt of outer edge
    dz00 = getD(rhos[0], h_grid)   # tilt of front outer edge
    dz0 = dz00                     # tilt of first leaf = tilt of outer edge
    
    for i in range(0,r-1): #z 
        # start with left edge of pixel, include right edge in loop
        # -> need to exclude last angle = right edge from looping   
        if (i+1)%(N_P[1]+1)==0: # if gap
            dz0 = getD(rhos[i+1], h_grid)  
            # do nothing
        else: 
            dz1 = getD(rhos[i+1], h_grid)  # tilt of right edge          
            for j in range(0,a-1): #x
                if (j+1)%(N_P[0]+1)==0: # if gap
                    dx0 = getD(alphas[j+1], h_grid)
                    # do nothing
                else:
                    dx1 = getD(alphas[j+1], h_grid)
                    L[index:index+4,:] = getBottomLines(X[j]+th_grid,X[j+1], 
                                                      Z[i]+th_grid,Z[i+1])
                    L[index+4:index+8,:] =getVerticalLines(X[j]+th_grid,X[j+1], 
                                                          Z[i]+th_grid,Z[i+1], 
                                                          dx0,dx1,dz0,dz1, 
                                                          h_grid)
                    L[index+8:index+12,:] = getTopLines(X[j]+th_grid,X[j+1], 
                                                      Z[i]+th_grid,Z[i+1], 
                                                      dx0,dx1,dz0,dz1,h_grid)
                    index = index + 12
                    dx0 = dx1 #tilt of le edge, reuse for ri edge of next pixel
            dz0 = dz1 # tilt of back edge, reuse for front edge of next pixel 
    
    #OUTER LINES:
    dx0 = dx00
    dz0 = dz00
    # dx1, dz1 can be directly reused    
    L[index:index+4,:] = getBottomLines (X[0],X[-1]+th_grid, 
                                         Z[0],Z[-1]+th_grid)
    L[index+4:index+8,:] = getVerticalLines(X[0],X[-1]+th_grid, Z[0],
                                            Z[-1]+th_grid, dx0,dx1,dz0,dz1, 
                                            h_grid)
    L[index+8:index+12,:] = getTopLines(X[0],X[-1]+th_grid, Z[0],Z[-1]+th_grid,
                                        dx0,dx1,dz0,dz1, h_grid)    
    return L'''




    


       
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        
