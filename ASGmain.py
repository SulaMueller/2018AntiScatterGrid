# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 17:49:47 2018

@author: Sula Mueller

DESCRIPTION:
    creates STLA or IGES file from parameters
    input: Input-Parameter-File.txt saved in same directory
    to create STLA or IGES file please do the following steps:
        I. INSERT values in Input_Parameter_File.txt
        II. SAVE file
        III. RUN this module
"""

import os
from smParam2File import readInputFromText, readAsMatrix, readAsIntMatrix
from smParam2File import createSTLAfromGridParams, createIGESfromGridParams

#-------------------------------------------------------------
''' I. READ INPUT FILE '''
thispath = os.path.dirname(os.path.realpath(__file__))
with open(thispath+'\Input-Parameter-File.txt','r') as inputfile:
    inputtext = inputfile.read()  # total file as string

#-------------------------------------------------------------
''' II. SET PARAMETERS '''

''' OUTPUT FILE '''
GRIDTYPE = readInputFromText(inputtext, 'GRIDTYPE')
if GRIDTYPE == '2D':
    grid = 0  # use grid as parameter specifying the gridtype
else:
    if (GRIDTYPE == 'x') or (GRIDTYPE == 'X'):
        grid = 1
    elif (GRIDTYPE == 'z') or (GRIDTYPE == 'Z'):
        grid = -1
    else:
        print('Error: enter dimension of 1D grid')
path = readInputFromText (inputtext, 'OUTpath')
if path == 'THIS':
    path = thispath
name = readInputFromText (inputtext, 'OUTname')
FILETYPE = readInputFromText (inputtext, 'FILETYPE')

''' GRID PARAMETERS '''
GR = float(readInputFromText (inputtext, 'GridRatio'))
th_grid = float(readInputFromText (inputtext, 'th_grid'))

''' SYSTEM PARAMETERS '''
SDD = float(readInputFromText (inputtext, 'SDD'))
EAD = readAsMatrix(inputtext, 'EdgeAxisDifference', grid)

''' MODULE '''
size_M = readAsMatrix(inputtext, 'size_M', 0)
                                         #  0: need matrix for either gridtype
edge_M0 = readAsMatrix(inputtext, 'edge_M0', grid)
edge_M1 = readAsMatrix(inputtext, 'edge_M1', grid)
gap_M = float(readInputFromText(inputtext, 'gap_M'))
n = int(readInputFromText(inputtext,'nnn'))

''' ELEMENT (TILE) '''
edge_E0 = readAsMatrix(inputtext, 'edge_T0', grid)
edge_E1 = readAsMatrix(inputtext, 'edge_T1', grid)
gap_E = readAsMatrix(inputtext, 'gap_T', grid)
N_E = readAsIntMatrix(inputtext, 'N_T', grid)

''' PIXEL '''
N_P = readAsIntMatrix(inputtext, 'N_P', grid)

if size_M == [0,0]:
    size_P = readAsMatrix(inputtext, 'size_P', grid)
    gap_P = readAsMatrix(inputtext, 'gap_P', grid)
    size_M = ([edge_M0[0] + edge_M1[0] + N_E[0]*(edge_E0[0] + edge_E1[0]
               + gap_E[0] + N_P[0]*(size_P[1]+gap_P[0])-gap_P[0])-gap_E[0],
               edge_M0[1] + edge_M1[1] + N_E[1]*(edge_E0[1] + edge_E1[1]
               + gap_E[1] + N_P[1]*(size_P[1]+gap_P[1]) -gap_P[1]) -gap_E[1]])
    print('size of M: ', str(size_M[0]), ' ', str(size_M[1]))

#-------------------------------------------------------------
''' III. CREATE OUTPUT FILE '''

if ((FILETYPE == 'STLA') or (FILETYPE == 'stla') or (FILETYPE == '.stl') or
    (FILETYPE == 'stl')):
    createSTLAfromGridParams(path, name, GR, th_grid, SDD,EAD, size_M,edge_M0,
                             edge_M1,gap_M, edge_E0,edge_E1,gap_E, N_E,N_P, n)

if ((FILETYPE == 'IGES') or (FILETYPE == 'iges') or (FILETYPE == '.igs') or
    (FILETYPE == 'igs')):
    createIGESfromGridParams(path, name, GR, th_grid, SDD, EAD, size_M,
                             edge_M0, edge_M1, gap_M, edge_E0, edge_E1, gap_E,
                             N_E, N_P, n)