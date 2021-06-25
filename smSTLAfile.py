# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 16:24:30 2018

DESCRIPTION: turns given input into STLA File
* input should be array created by smFacets
    -> array with normals and vertice saved as linear array

@author: Sula Mueller
"""

from decimal import Decimal

''' turns float into string (exp format, as requested by stla file reader)
REMARK: '%.5E' specifies number of digits/ precision '''
def formatIntoDecimalString(number):
    return str('%.4E'%Decimal(number))

''' createSTLAFile
DESCRIPTION: turns given input into STLA File
INPUT:
    * 1D array with all values stored in linear fashion
        -> array with normals and vertice saved as linear array
    * path: location to save created file (eg: 'C:/ASG/gridfiles/')
        end with / or \
    * name: of file without extension, eg: 'example'
OUTPUT: .stl file saved at specified location '''
def createSTLAFile(path, name, Facets):
    path = path + '\\'
    file = open(path + name + '.stl', 'w')  # initially open file
    ''' HEAD: '''
    file.write('solid ' + name +'\n')
    ''' CONTENT: '''
    f = len(Facets)
    print(str(int(f/12)), 'facets')  # to inform user ^^
    ''' FORMAT ALL THE NUMBERS: '''
    F = []
    for i in range(0, f):
        F.append(formatIntoDecimalString(Facets[i]))
    ''' WRITE IN FILE: '''
    for i in range(0, int(f/12)):  # 12 entries for 1 facet
        file.write('  facet normal '+F.pop(0)+' '+F.pop(0)+' '+F.pop(0) +'\n'
                   +'    outer loop'+'\n'+'      vertex '+ F.pop(0)+' '
                   +F.pop(0)+' '+F.pop(0) +'\n'+'      vertex '+ F.pop(0)+' '
                   +F.pop(0)+' '+F.pop(0) +'\n'+'      vertex '+ F.pop(0)+' '
                   +F.pop(0)+' '+F.pop(0)+'\n'+'    endloop'+'\n'+'  endfacet'
                   +'\n')
    ''' FOOT: '''
    file.write('endsolid')
    file.close()
    print('Done')
    return