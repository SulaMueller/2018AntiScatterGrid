# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 07:12:22 2018

DESCRIPTION: turns given input into IGES File
INPUT: Angles as determined by smAngles -> determineAllAnglesAndCoordinates
OUTPUT: .igs file

@author: Sula Mueller
"""

import time
from smPsection import allP_Entries

''' creates 7 digit string from value x and enough spaces (right justified) '''
# eigth in block is section entry ('D', or 'P')
def spaced7(x):
    xx = str(x)
    return (7-len(xx))*' ' + xx
 
''' creates 8 digit string from value x and enough spaces (right justified) '''
def spaced8(x):
    xx = str(x)
    return (8-len(xx))*' ' + xx

''' creates string length ls from value x and enough spaces 
    (right justified) '''
def spacedx(x, ls):
    xx = str(x)
    return (ls-len(xx))*' ' + xx

def getDateAndTime():  # needed for G section
    t = time.localtime(time.time())
    
    if t.tm_mon<10:
        m = '0' + str(t.tm_mon)
    else:
        m = str(t.tm_mon)
    
    if t.tm_mday<10:
        d = '0' + str(t.tm_mday)
    else:
        d = str(t.tm_mday)
    
    if t.tm_hour<10:
        h = '0' + str(t.tm_hour)
    else:
        h = str(t.tm_hour)
    
    if t.tm_min<10:
        mn = '0' + str(t.tm_min)
    else:
        mn = str(t.tm_min)
    
    if t.tm_sec<10:
        s ='0' + str(t.tm_sec)
    else:
        s = str(t.tm_sec)
    return str(t.tm_year) + m + d + '.' + h + mn + s

''' createIGESFile
DESCRIPTION: creates IGES file name.igs, saves at path, as specified by .igs 
             guidelines
INPUT:
    Angles (from smAngles.determineAllAnglesAndCoordinates), Grid Ratio,
    thickness of grid leaves
OUTPUT: file '''
def createIGESFile(path, name, Angles, GR, th_grid):
    file = open(path + '/' + name + '.igs', 'w')
    ''' S SECTION: '''
    file.write(72*' ' + 'S      1' + '\n')
    
    ''' G SECTION: '''
    n = len(name) + 4  # +4 because .igs
    NAME = name.upper()
    p = len(path) + n
    dt = getDateAndTime()
    
    g = ('1H,,1H;,'+str(n)+'H'+NAME+'.IGS,'+str(p)+'H'+path+name+
         '.igs,10HASGCreator,3H1.1,32,8,23,11,52,'+str(n)+'H'+NAME+
         '.IGS,1.,1,2HMM,4,0.7,15H'+dt+
         ',0.1000000000E-003,4747.868164,,,11,,,;')
    g = g + (72-len(g)%72)*' '  # to fill up last line with spaces
    len_g = int(len(g)/72)  # number of lines g fills
    for i in range(0, len_g):  # each line
        file.write(g[i*72:i*72+71] + ' G      ' + str(i+1) + '\n')
    
    ''' CREATE P SECTION: '''
    # needs to be created first because pointers point to entities
    # (but write to file later; after D section)
    # use functions defined in smPsection
    P = allP_Entries(Angles, GR, th_grid)  # from smPsection.py
        # P[0] -> first nextline = 1
        # P[1+3i] -> string as written in P section
        # P[2+3i] -> type as integer
        # P[3+3i] -> line index for NEXT entry ("nextline")
    lenP = int((len(P)-1)/3)  # number of actual P entries
    # -1 since first entry initalized as 1 (= first "nextline" for first line)
    # /3 since alternating strings, "nextlines", types
    print ('\n' + 'number of P entries: ' + str(lenP))
    
    ''' D SECTION: '''
    p = spaced8(1) # pointer for entry in P section (P[1] itnitialized as '1')
    for i in range(0, lenP):
        linesp = spaced8(P[3*i+3]-P[3*i])
                 # number of lines in P for this entity
        x = spaced7(2*i+1)  # line index of current line
        y = spaced7(2*i+2)  # line index of 2nd line (its a 2 line paragraph)
        typ = spaced8(P[3*i+2])
        file.write(typ + p +
                   '       0       1       0       0       0       0       0D'
                   + x + '\n')
        file.write(typ + '       1' + '       4' + linesp +
                   '       0                               1D' + y + '\n')
        p = spaced8(P[3*i+3])  # get pointer for next line as nextline
    
    ''' WRITE P SECTION: '''
    for i in range(0, lenP):
        file.write(P[3*i+1])
    
    ''' T SECTION: '''
    file.write('S      1G      ' + str(len_g) + 'D' + spaced7(lenP*2) + 'P'
               + spaced7(P[-1]-1) + 40*' ' + 'T      1')
    file.close()
    print('Done')
    return