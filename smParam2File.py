# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 14:55:30 2018

@author: Sula Mueller

DESCRIPTION:
    defines functions used by ASGmain
"""

from smAngles import determineAnglesSparseInfo
from smAngles import determineAllAnglesAndCoordinates
from smFacets import FacetizeEntireGrid
from smSTLAfile import createSTLAFile
from smIGESfile import createIGESFile

''' readInputFromText
DESCRIPTION: extracts value/ input from given file-string by designation of
    value
INPUT:
    * inputtext: content of Input_Parameter_File as string
    * valuename: name of value that gets extracted (eg 'FILETYPE' or 'SDD')
OUTPUT:
    * value of designated valuename (eg 'STLA' or '120')
    * output always as string '''
def readInputFromText(inputtext, valuename):
    i = inputtext.find(valuename)  # index of valuename
    i = inputtext.find('=', i)   # index of = after valuename
    substring = inputtext[i:-1]
    value = substring.split()  # returns array of all non-space entries
    return value[1]  # first entry is '=', second entry should be wanted value

''' readAsMatrix
DESCRIPTION: extracts value/ input from given file-string by designation of
    value
INPUT:
    * inputtext: content of Input_Parameter_File as string
    * valuename: name of value that gets extracted (eg 'N_P')
    * OneDgrid: 0 -> 2D grid
                1 -> 1D with leaves along x-axis
                -1 -> 1D with leaves along z-axis
OUTPUT:
    * value of designated valuename (eg [32,32])
    * output always as 1x2 matrix type float
    * if 1D grid: values of OTHER dimension are set to 0 '''
def readAsMatrix(inputtext, valuename, gridtype):
    if gridtype == 1:  # 1D with leaves along x-axis
        val = float(readInputFromText(inputtext, valuename))
        return [0,val]
    if gridtype == -1:  # 1D with leaves along z-axis
        val = float(readInputFromText(inputtext, valuename))
        return [val,0]
    # else: 2D grid
    i = inputtext.find(valuename)  # index of valuename
    i = inputtext.find('[', i)  # index of first [ after valuename
    j = inputtext.find(']', i)
    substring = inputtext[i+1:j]  # extract string between [] (without [])
    if substring == '[]':
        return ('[0,0]')
    substring.replace(';', ',')
    value = substring.split(',')
    valu = value[0]  # first entry
    valu.strip()  # remove ' '
    x = float(valu)
    valu = value[1]  # second entry
    valu.strip()  # remove ' '
    z = float(valu)
    return [x,z]

''' readAsIntMatrix
same as readAsMatrix but output is matrix of integers
if 1D grid: values of OTHER dimension are set to 1
(readAsIntMatrix gets only called by N_E, N_P -> treat as if 1 pixel) '''
def readAsIntMatrix(inputtext, valuename, gridtype):
    if gridtype == 1:  # 1D with leaves along x-axis
        val = int(readInputFromText(inputtext, valuename))
        return [1,val]
    if gridtype == -1:  # 1D with leaves along z-axis
        val = int(readInputFromText(inputtext, valuename))
        return [val,1]
    # else: 2D grid
    i = inputtext.find(valuename)  # index of valuename
    i = inputtext.find('[', i)  # index of first [ after valuename
    j = inputtext.find(']', i)
    substring = inputtext[i+1:j]  # extract string between [] (without [])
    if substring == '[]':
        return ('[0,0]')
    substring.replace(';',',')
    value = substring.split(',')
    valu = value[0]  # first entry
    valu.strip()  # remove ' '
    x = int(valu)
    valu = value[1]  # second entry
    valu.strip()  # remove ' '
    z = int(valu)
    return [x,z]

''' does same as createSTLAfromGridParams with simpler geometry '''
def createSTLAfromGridParams_Sparse(path, name, GR, th_grid, SDD, delta_M, N_E,
                                    N_P, n):
    Angles = determineAnglesSparseInfo(SDD, delta_M, N_E, N_P, n)
    Facets = FacetizeEntireGrid(Angles, GR, th_grid)
    createSTLAFile(path, name, Facets)
    return

def createSTLAfromGridParams(path, name, GR, th_grid, SDD, EAD, size_M,
                             edge_M0, edge_M1, gap_M, edge_E0, edge_E1, gap_E,
                             N_E, N_P, n):
    Angles = determineAllAnglesAndCoordinates(SDD, EAD, size_M, edge_M0,
                                              edge_M1, gap_M, edge_E0,
                                              edge_E1, gap_E, N_E, N_P, n)
    Facets = FacetizeEntireGrid(Angles, GR, th_grid)
    createSTLAFile(path, name, Facets)
    return
    
def createIGESfromGridParams(path, name, GR, th_grid, SDD, EAD, size_M,
                             edge_M0, edge_M1, gap_M, edge_E0, edge_E1, gap_E,
                             N_E, N_P, n):
    Angles = determineAllAnglesAndCoordinates(SDD, EAD, size_M, edge_M0,
                                              edge_M1, gap_M, edge_E0, edge_E1,
                                              gap_E, N_E, N_P, n)
    createIGESFile(path, name, Angles, GR, th_grid)
    return