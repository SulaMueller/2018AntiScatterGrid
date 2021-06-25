# -*- coding: utf-8 -*-

"""
Created on Tue May  8 17:13:01 2018

DESCRIPTION: creates P section for IGES File
INPUT: Angles as determined by smAngles -> determineAllAnglesAndCoordinates
MAIN ROUTINE: allP_Entries(Angles, GR, th_grid)

REMARK: nextline, entry_index
from specifications of IGES file:
* entry_index:
    * specifies the entry, is identical for entire P entry
    * independent of how many lines have been used up for that entry
    * only uses odd numbers (hence +2 for each entry)
    * entry_index can be used as "pointer" (DE) to that entry
* line counter ("nextline"):
    * continuously counts used up lines
    * independent of which P entry it is used for
    * set variable "nextline" = number of used up lines
    -> "nextline" is line index for next entry (needed in D-section)

SECTIONS:
    0: Basics (P_Entry, Points, Lines)
    I: Surfaces directly boardering Pixels
    II: Edge faces of entire grid
    III: Bottom and Top surface
        A: between pixels ("alpha")
        B: between pixel rows ("rho")
    IV: Main routine

@author: Sula Mueller
"""

'''---------------------------------------------------------------------------
SECTION 0: BASICS '''

import numpy as np
from smFacets import getD

''' creates 7 digit string from value x and enough spaces (right justified) '''
# eigth in block is section entry ('D', or 'P')
def spaced7(x):
    xx = str(x)
    return (7-len(xx))*' ' + xx

''' creates 8 digit string from value x and enough spaces (right justified) '''
def spaced8(x):
    xx = str(x)
    return (8-len(xx))*' ' + xx

''' reverses direction of line '''
def reverseLine(L):
    return [L[3], L[4], L[5], L[0], L[1], L[2]]

''' P_Entry
DESCRIPTION: creates one "entry" for P section from given parameters
INPUT: entry_index -> index of entry (independent of used up lines)
       nextline -> line counter for first line
       typ -> type of entity as integer
       params -> parameters of entity (if point: [x,y,z], if line: P1, P2)
       pointFlag -> true, if point, line or curve, false if surface
OUTPUT:
    P, alternating P[0+3i] -> string as written in P section (incl. '/n')
                   P[1+3i] -> type as integer
                   P[2+3i] -> line index for NEXT entry ("nextline") '''
def P_Entry(entry_index, nextline, typ, params, pointFlag):
    nextline = int(nextline)
    p = str(typ) + ','  # p: string output without '\n'
    for i in range(0, len(params[:])):
        p = p + str(params[i]) + ','  # float with 4 digits
    if pointFlag:
        p = p + '0,0;'  # line scheme for points, lines and curves
    else:
        p = p + '0,1,43;'  # line scheme for surfaces
    string = ''  # string output with added '\n'
    estring = spaced8(entry_index) + 'P'  # string form of entry_index
                                          # gets added at each line
    if len(p)>64:
        while (len(p)>64):
            # slice after next comma before 64 AND REMOVE:
            ind = p.rfind(',', 0, 64)
            pp = p[0:ind+1]  # sliced off part
            p = p[ind+1:len(p)]  # leftover part
            string = (string + pp + (64-len(pp))*' ' + estring
                      + spaced7(nextline) + '\n')
            nextline = int(nextline + 1)
    # last part (or if len < 64 in the first place):
    string = string + p + (64-len(p))*' ' + estring + spaced7(nextline) + '\n'
    return [string, typ, int(nextline+1)]

''' RuledSurface
DESCRIPTION: returns surface as P_Entry
INPUT: L -> pointers of 2 antiparallel lines on that surface
OUTPUT: 1 surface as P_Entry '''
def RuledSurface(L, entry_index, nextline):
    params = [L[0], L[1], 1, 0]
    pointFlag = False
    return P_Entry(entry_index, nextline, 118, params, pointFlag)

''' BoundedSurfaces
DESCRIPTION: gets P_Entries of bounded surfaces
INPUT: SP -> pointers to unbounded surfaces
       BP -> pointers to boundaries
OUTPUT: P_Entries '''
def BoundedSurfaces(SP, BP, entry_index, nextline):
    pointFlag = False
    lenP = len(SP)  # assume, len(SP) = len(BP)
    P = []  # output
    for i in range(0,lenP):
        params = [0,SP[i],1,BP[i]]
        P = P + P_Entry(entry_index, nextline, 143, params, pointFlag)
        nextline = int(P[-1])
        entry_index = entry_index + 2  # per definition of P section
    return P

''' BoundedSurface
DESCRIPTION: gets P_Entry of 1 bounded surface
INPUT: SP -> pointer to unbounded surface
       BP -> pointer to boundary
OUTPUT: P_Entry '''
def BoundedSurface(SP, BP, entry_index, nextline):
    pointFlag = False
    params = [0,SP,1,BP]
    return P_Entry(entry_index, nextline, 143, params, pointFlag)

'''---------------------------------------------------------------------------
SECTION I: PIXELS '''

''' PixelPoints
DESCRIPTION: returns line array of all points surrounding pixels
INPUT: coordinates of 4 bottom points ([x0,z0], [x0,z1], [x1,z1], [x1,z0])
           (x0, z0 should be +th_grid)
       dx, dz -> tilt of upper points
       h_grid -> height of grid
OUTPUT: line array [8,3] '''
def PixelPoints(x0, x1, z0, z1, dx0, dx1, dz0, dz1, h_grid):
    P = np.zeros((8,3))
    P[0,:] = ["%.4f"%x0, 0, "%.4f"%z0]  # "%.4f"% to limit to 4 digits
    P[1,:] = ["%.4f"%x1, 0, "%.4f"%z0]
    P[2,:] = ["%.4f"%x1, 0, "%.4f"%z1]
    P[3,:] = ["%.4f"%x0, 0, "%.4f"%z1]
    P[4,:] = ["%.4f"%(x0+dx0), h_grid, "%.4f"%(z0+dz0)]
    P[5,:] = ["%.4f"%(x1+dx1), h_grid, "%.4f"%(z0+dz0)]
    P[6,:] = ["%.4f"%(x1+dx1), h_grid, "%.4f"%(z1+dz1)]
    P[7,:] = ["%.4f"%(x0+dx0), h_grid, "%.4f"%(z1+dz1)]
    return P

''' PixelLines
DESCRIPTION: returns line array of all lines surrounding pixels
INPUT: coordinates of 4 bottom points ([x0,z0], [x0,z1], [x1,z1], [x1,z0])
           (x0, z0 should be +th_grid)
       dx, dz -> tilt of upper points
       h_grid -> height of grid
OUTPUT: line array [12,6]
REMARKS:
    * DON'T do for last alpha,rho since need i+1
    * "%.4f"% to round on 4 digits '''
def PixelLines(x0, x1, z0, z1, dx0, dx1, dz0, dz1, h_grid):
    L = np.zeros((12,6))
    ''' BOTTOM: '''
    L[0,:] = ["%.4f"%x0,0,"%.4f"%z0, "%.4f"%x1,0,"%.4f"%z0]
    L[1,:] = ["%.4f"%x1,0,"%.4f"%z0, "%.4f"%x1,0,"%.4f"%z1]
    L[2,:] = ["%.4f"%x1,0,"%.4f"%z1, "%.4f"%x0,0,"%.4f"%z1]
    L[3,:] = ["%.4f"%x0,0,"%.4f"%z1, "%.4f"%x0,0,"%.4f"%z0]
    ''' VERTICAL: '''
    L[4,:] = ["%.4f"%x0,0,"%.4f"%z0, "%.4f"%(x0+dx0),h_grid,"%.4f"%(z0+dz0)]
    L[5,:] = ["%.4f"%x1,0,"%.4f"%z0, "%.4f"%(x1+dx1),h_grid,"%.4f"%(z0+dz0)]
    L[6,:] = ["%.4f"%x1,0,"%.4f"%z1, "%.4f"%(x1+dx1),h_grid,"%.4f"%(z1+dz1)]
    L[7,:] = ["%.4f"%x0,0,"%.4f"%z1, "%.4f"%(x0+dx0),h_grid,"%.4f"%(z1+dz1)]
    ''' TOP: '''
    L[8,:] = ["%.4f"%(x0+dx0),h_grid,"%.4f"%(z0+dz0), "%.4f"%(x1+dx1),h_grid,
              "%.4f"%(z0+dz0)]
    L[9,:] = ["%.4f"%(x1+dx1),h_grid,"%.4f"%(z0+dz0), "%.4f"%(x1+dx1),h_grid,
              "%.4f"%(z1+dz1)]
    L[10,:] = ["%.4f"%(x1+dx1),h_grid,"%.4f"%(z1+dz1), "%.4f"%(x0+dx0),h_grid,
               "%.4f"%(z1+dz1)]
    L[11,:] = ["%.4f"%(x0+dx0),h_grid,"%.4f"%(z1+dz1), "%.4f"%(x0+dx0),h_grid,
               "%.4f"%(z0+dz0)]
    return L

''' PixelRuledSurfaces
DESCRIPTION: returns array of all surfaces boardering a pixel
INPUT: LP -> pointers to all lines boardering pixel
OUTPUT: 4 surfaces as P_Entries '''
def PixelRuledSurfaces(LP, entry_index, nextline):
    pointFlag = False
    # front:
    params = [LP[0],LP[8],0,0]
    P = P_Entry(entry_index, nextline, 118, params, pointFlag)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # right:
    params = [LP[1],LP[9],0,0]
    P = P + P_Entry(entry_index, nextline, 118, params, pointFlag)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # back:
    params = [LP[2],LP[10],0,0]
    P = P + P_Entry(entry_index, nextline, 118, params, pointFlag)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # left:
    params = [LP[3],LP[11],0,0]
    return P + P_Entry(entry_index, nextline, 118, params, pointFlag)

''' PixelBoundaries
DESCRIPTION: gets the boundaries of 4 pixel boardering surfaces from Lines
INPUT: LP -> pointers to all lines boardering pixel
       SP -> pointers to 4 unbounded surfaces
OUTPUT: 4 boundaries for 1 pixel as P_Entries '''
def PixelBoundaries(LP, SP, entry_index, nextline):
    pointFlag = True
    # front:
    params = [0,0,SP[0],4,LP[0],0,0,LP[4],0,0,LP[8],1,0,LP[7],1,0]
    P = P_Entry(entry_index, nextline, 141, params, pointFlag)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # right:
    params = [0,0,SP[1],4,LP[1],0,0,LP[5],0,0,LP[9],1,0,LP[4],1,0]
    P = P + P_Entry(entry_index, nextline, 141, params, pointFlag)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # back:
    params = [0,0,SP[2],4,LP[2],0,0,LP[6],0,0,LP[10],1,0,LP[5],1,0]
    P = P + P_Entry(entry_index, nextline, 141, params, pointFlag)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # left:
    params = [0,0,SP[3],4,LP[3],0,0,LP[7],0,0,LP[11],1,0,LP[6],1,0]
    return P + P_Entry(entry_index, nextline, 141, params, pointFlag)

''' P_EntriesOnePixel
DESCRIPTION:
    creates ALL P_Entries of points, lines and bounded surfaces boardering one
    pixel
INPUT:
    values -> x0, x1, z0, z1, dx0, dx1, dz0, dz1, h_grid
    = X[j]+th_grid, X[j+1], Z[i]+th_grid, Z[i+1], dx0, dx1, dz0, dz1, h_grid
OUTPUT: 32 P_Entries:
    4 bottom points
    4 top points
    4 bottom lines
    4 vertical lines
    4 top lines
    4 ruled, unbounded surfaces (created by sweeping bottom, top lines)
    4 boundaries of those (created by lines)
    4 bounded surfaces (created by surfaces and corresponding boundaries)
1 P_Entry = 3 array elements:
    P[0+3i] -> string as written in P section
    P[1+3i] -> type as integer
    P[2+3i] -> line index for NEXT entry ("nextline")
LP -> pointers to all lines boardering pixel '''
def P_EntriesOnePixel(entry_index, nextline, values):
    P = []  # output
    LP = np.zeros(12, dtype=int)  # pointers for lines from entry_index
    ''' points as [8,3] '''
    Points = PixelPoints(values[0],values[1],values[2],values[3],values[4],
                         values[5],values[6],values[7],values[8])
    pointFlag = True
    for i in range(0,8):  # write Points as P_Entries:
        P = P + P_Entry(entry_index, nextline, 116, Points[i], pointFlag)
        nextline = int(P[-1])
        entry_index = entry_index + 2  # per definition of P section
    ''' lines as [12,6] '''
    Lines = PixelLines(values[0],values[1],values[2],values[3],values[4],
                       values[5],values[6],values[7],values[8])
    for i in range(0,12):  # write Lines as P_Entries:
        LP[i] = entry_index
        P = P + P_Entry(entry_index, nextline, 110, Lines[i], pointFlag)
        nextline = int(P[-1])
        entry_index = entry_index + 2  # per definition of P section
    ''' ruled surfaces '''
    P = P + PixelRuledSurfaces(LP, entry_index, nextline)
    e = entry_index
    SP = [e, e+2, e+4, e+6]  # pointers for surfaces as entry_indice
    nextline = int(P[-1])
    entry_index = entry_index + 8  # +2*(4 surfaces)
    ''' boundaries for inner surfaces '''
    #P = P + PixelBoundaries(LP, SP, entry_index, nextline)
    e = entry_index
    BP = [e, e+2, e+4, e+6]  # pointers for boundaries as entry_indice
    nextline = int(P[-1])
    entry_index = entry_index + 8  # +2*(4 boundaries)
    ''' bounded surfaces '''
    P = P + BoundedSurfaces(SP, BP, entry_index, nextline)
    return [P, LP]

'''---------------------------------------------------------------------------
SECTION II: EDGES OF GRID '''

''' OuterEdgesLines
DESCRIPTION:
    creates P_Entries of points and front/back lines of outer edges of grid
    equivalent to P_EntriesOnePixel, but exclude left, right lines
        (otherwise would have redundant definition of lines)
        (boundaries & surfaces defined in OuterEdgesSurfaces)
INPUT:
    values -> x0, x1, z0, z1, dx0, dx1, dz0, dz1, h_grid
OUTPUT:
    P -> P_entries of points and lines, except left/ right lines
    oeL ->  pointers to [bottom front+back, vertical, top front+back] lines '''
def OuterEdgesLines(entry_index, nextline, values):
    P = []  # output
    oeL = np.zeros(8, dtype=int)  # pointers for lines from entry_index
    ''' points as [8,3] '''
    Points = PixelPoints(values[0],values[1],values[2],values[3],values[4],
                         values[5],values[6],values[7],values[8])
    pointFlag = True
    for i in range(0,8):  # write Points as P_Entries:
        P = P + P_Entry(entry_index, nextline, 116, Points[i], pointFlag)
        nextline = int(P[-1])
        entry_index = entry_index + 2  # per definition of P section
    ''' lines as [8,6] '''
    Lines = PixelLines(values[0],values[1],values[2],values[3],values[4],
                       values[5],values[6],values[7],values[8])
    Lines = [Lines[0], Lines[2], Lines[4], Lines[5], Lines[6], Lines[7],
             Lines[8], Lines[10]]  # exclude left/ right edges
    for i in range(0,8):  # write Lines as P_Entries:
        oeL[i] = entry_index
        P = P + P_Entry(entry_index, nextline, 110, Lines[i], pointFlag)
        nextline = int(P[-1])
        entry_index = entry_index + 2  # per definition of P section
    return [P, oeL]

''' OuterEdgesRuledSurf
DESCRIPTION: returns array of outer edge surfaces of grid
INPUT: oeL -> pointers to all lines [bottom front+back, vert, top front+back]
OUTPUT: 4 surfaces as P_Entries '''
def OuterEdgesRuledSurfs(oeL, entry_index, nextline):
    pointFlag = False
    # front:
    params = [oeL[2],oeL[3],0,0]
    P = P_Entry(entry_index, nextline, 118, params, pointFlag)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # right:
    params = [oeL[3],oeL[4],0,0]
    P = P + P_Entry(entry_index, nextline, 118, params, pointFlag)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # back:
    params = [oeL[4],oeL[5],0,0]
    P = P + P_Entry(entry_index, nextline, 118, params, pointFlag)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # left:
    params = [oeL[5],oeL[2],0,0]
    return P + P_Entry(entry_index, nextline, 118, params, pointFlag)

''' OuterEdgesRLboundary
DESCRIPTION: gets boundary for right or left surface of outer edges of grid
INPUT: vLP -> pointers to vertical lines (bottom to top)
       bLP -> pointers to all bottom lines (many partials, front to back)
       tLP -> pointers to all top lines (many partials, front to back)
       SP -> pointers to unbounded surfaces
OUTPUT: boundary for outer edge face as P_Entry
REMARK: can also be used for [bottom, top] surfaces between pixel rows '''
def OuterEdgesRLboundary(vLP, bLP, tLP, SP, entry_index, nextline):
    pointFlag = True
    lenC = 2  # number of partial curves, initialize with 2 vertical lines
    ''' vertical 1: '''
    params =  [vLP[0],1,0]
    ''' bottom: '''
    if isinstance(bLP, np.int32):
        params = np.concatenate((params, [bLP,0,0]))
        lenC = lenC + 1
    else:
        lenLP = len(bLP)
        lenC = lenC + lenLP
        for i in range(0, lenLP):
            params = np.concatenate((params, [bLP[i],0,0]))
    ''' vertical 2: '''
    params = np.concatenate((params, [vLP[1],0,0]))
    ''' top: '''
    if isinstance(tLP, np.int32):
        params = np.concatenate((params, [tLP,0,0]))
        lenC = lenC + 1
    else:
        lenLP = len(tLP)
        lenC = lenC + lenLP
        for i in range(0, lenLP):
            params = np.concatenate((params, [tLP[-i],1,0]))
    ''' first part of boundary definition: '''
    params = np.concatenate(([0,0,SP,lenC], params))
    return P_Entry(entry_index, nextline, 141, params, pointFlag)

''' OuterEdgesBoundaries
DESCRIPTION: gets the boundaries of 4 outer edge surfaces of grid
INPUT: oeL -> pointers to all lines except left/ right lines
       l/r b/t L,  -> pointers to left/right bottom/top lines of edges
       SP -> pointers to 4 unbounded surfaces
OUTPUT: 4 boundaries for outer edges as P_Entries '''
def OuterEdgesBoundaries(oeL, lbL, rbL, ltL, rtL, SP, entry_index, nextline):
    pointFlag = True
    # front:
    params = [0,0,SP[0],4,oeL[0],0,0,oeL[2],0,0,oeL[6],1,0,oeL[5],1,0]
    P = P_Entry(entry_index, nextline, 141, params, pointFlag)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # right:
    vL = [oeL[2], oeL[3]]  # pointers to vertical lines of right edge
    P = P + OuterEdgesRLboundary(vL, rbL, rtL, SP[1], entry_index, nextline)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # back:
    params = [0,0,SP[2],4,oeL[1],0,0,oeL[4],0,0,oeL[7],1,0,oeL[3],1,0]
    P = P + P_Entry(entry_index, nextline, 141, params, pointFlag)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # left:
    vL = [oeL[4], oeL[5]]  # pointers to vertical lines of left edge
    return (P + OuterEdgesRLboundary(vL, lbL, ltL, SP[3], entry_index,
                                     nextline))

''' OuterEdgesSurfaces
DESCRIPTION: wraps up surface routines to create 4 outer edge surfaces of grid
INPUT: oeLP -> pointers of all lines except left/ right lines
       oelrLP -> pointers to left/ rigth lines of edges
OUTPUT: four bounded surfaces of outer edges as P_Entries '''
def OuterEdgesSurfaces(oeL, edgeRL, SP, entry_index, nextline):
    lenE = len(edgeRL)
    lbL = edgeRL[0:lenE:4]  # pointers to left bottom partial lines
    rbL = edgeRL[1:lenE:4]  # pointers to right bottom partial lines
    ltL = edgeRL[2:lenE:4]  # pointers to left top partial lines
    rtL = edgeRL[3:lenE:4]  # pointers to right top partial lines
    ''' ruled surfaces '''
    P = OuterEdgesRuledSurfs(oeL, entry_index, nextline)
    e = entry_index
    SP = [e, e+2, e+4, e+6]  # pointers for surfaces
    nextline = int(P[-1])
    entry_index = entry_index + 8  # +2*(4 surfaces)
    ''' boundaries for surfaces '''
    P = P + OuterEdgesBoundaries(oeL, lbL, rbL, ltL, rtL, SP, entry_index,
                                 nextline)
    e = entry_index
    BP = [e, e+2, e+4, e+6]  # pointers for boundaries
    nextline = int(P[-1])
    entry_index = entry_index + 8  # +2*(4 boundaries)
    ''' bounded surfaces '''
    return P + BoundedSurfaces(SP, BP, entry_index, nextline)

'''---------------------------------------------------------------------------
SECTION III: BOTTOM & TOP '''

''' PART A: ALPHA '''

''' BotTopAlpPartialLineParams
DESCRIPTION: gets additional lines surrounding [bottom, top] surfaces between
             pixels
INPUT: coordinates of 4 bottom points ([x0,z0], [x0,z1], [x1,z1], [x1,z0]):
       [X[j], X[j]+th_grid, Z[i]+th_grid, Z[i+1]]
       dx, dz -> tilt of upper points
       h_grid -> height of grid
OUTPUT: line array [4,6] '''
def BotTopAlpPartialLineParams(x0, x1, z0, z1, dx0, dx1, dz0, dz1, h_grid):
    L = np.zeros((4,6))
    L[0,:] = ["%.4f"%x0,0,"%.4f"%z0, "%.4f"%x1,0,"%.4f"%z0]
    L[1,:] = ["%.4f"%x0,0,"%.4f"%z1, "%.4f"%x1,0,"%.4f"%z1]
    L[2,:] = ["%.4f"%(x0+dx0),h_grid,"%.4f"%(z0+dz0), "%.4f"%(x1+dx1),h_grid,
              "%.4f"%(z0+dz0)]
    L[3,:] = ["%.4f"%(x0+dx0),h_grid,"%.4f"%(z1+dz1), "%.4f"%(x1+dx1),h_grid,
              "%.4f"%(z1+dz1)]
    return L

''' BotTopAlpBoundaries
DESCRIPTION: gets the boundaries of [bottom, top] surfaces between pixels
INPUT: parL -> pointers to partial lines between pixels
       alps -> pointers to [left bottom, right bottom, l top, r top] lines of
               space in between
       SP -> pointers to unbounded surfaces
OUTPUT: 2 boundaries as P_Entries '''
def BotTopAlpBoundaries(parL, alps, SP, entry_index, nextline):
    pointFlag = True
    # bottom:
    params = [0,0,SP[0],4,parL[0],1,0,alps[0],0,0,parL[1],0,0,alps[1],0,0]
    P = P_Entry(entry_index, nextline, 141, params, pointFlag)
    nextline = int(P[-1])
    entry_index = entry_index + 2  # per definition of P section
    # top:
    params = [0,0,SP[1],4,parL[2],1,0,alps[2],0,0,parL[3],0,0,alps[3],0,0]
    return P + P_Entry(entry_index, nextline, 141, params, pointFlag)

''' BotTopAlphaSurfaces
DESCRIPTION: gets all P_entries of [bottom, top] surfaces between pixels
INPUT: alps -> pointers to [left bottom, right bottom, l top, r top] lines of
               space in between
       SP -> pointers to unbounded surfaces
       values -> [X[j], X[j]+th_grid, Z[i]+th_grid, Z[i+1], dx0, dx1, dz0, dz1,
                  h_grid]
OUTPUT: P -> all P_Entries for both surfaces
        parLP -> pointers to additional partial lines between pixels '''
def BotTopAlphaSurfaces(alps, SP, entry_index, nextline, values):
    P = []  # output
    parLP = np.zeros(4, dtype=int)  # pointers for lines from entry_index
    ''' additional partial lines: '''
    parL = BotTopAlpPartialLineParams(values[0],values[1],values[2],values[3],
                                      values[4],values[5],values[6],values[7],
                                      values[8])
    pointFlag = True
    for i in range(0,4):  # write parL as P_Entries:
        parLP[i] = entry_index
        P = P + P_Entry(entry_index, nextline, 110, parL[i], pointFlag)
        nextline = int(P[-1])
        entry_index = entry_index + 2  # per definition of P section
    ''' boundaries: '''
    BP = [entry_index, entry_index+2]  # pointers for boundaries
    P = P + BotTopAlpBoundaries(parLP, alps, SP, entry_index, nextline)
    nextline = int(P[-1])
    entry_index = entry_index + 4  # + 2*(2 boundaries)
    ''' bounded surfaces '''
    P = P + BoundedSurfaces(SP, BP, entry_index, nextline)
    return [P, parLP]

''' PART B: RHO '''

''' BotTopRhoPartialPointParams
DESCRIPTION: returns line array of 4 additional points for partial lines
INPUT: coordinates of 2 bottom points ([x0,z1], [x1,z1])
       [X[0], X[-1]+th_grid, Z[i]+th_grid]
       dx, dz -> tilt of upper points
       h_grid -> height of grid
OUTPUT: line array [4,3] '''
def BotTopRhoPartialPointParams(x0, x1, z1, dx0, dx1, dz1, h_grid):
    P = np.zeros((4,3))
    P[0,:] = ["%.4f"%x1, 0, "%.4f"%z1]
    P[1,:] = ["%.4f"%x0, 0, "%.4f"%z1]
    P[2,:] = ["%.4f"%(x1+dx1), h_grid, "%.4f"%(z1+dz1)]
    P[3,:] = ["%.4f"%(x0+dx0), h_grid, "%.4f"%(z1+dz1)]
    return P

''' BotTopRhoPartialLineParams
DESCRIPTION: gets [bl, br, tl, tr] edge lines of surfaces between pixel rows
INPUT: coordinates of 4 bottom points ([x0,z0], [x0,z1], [x1,z1], [x1,z0]):
       [X[0], X[-1]+th_grid, Z[i], Z[i]+th_grid]
       dx, dz -> tilt of upper points
       h_grid -> height of grid
OUTPUT: line array [4,6] '''
def BotTopRhoPartialLineParams(x0, x1, z0, z1, dx00, dx11, dz0, dz1, h_grid):
    L = np.zeros((4,6))
    L[0,:] = ["%.4f"%x0,0,"%.4f"%z0, "%.4f"%x0,0,"%.4f"%z1]
    L[1,:] = ["%.4f"%x1,0,"%.4f"%z0, "%.4f"%x1,0,"%.4f"%z1]
    L[2,:] = ["%.4f"%(x0+dx00),h_grid,"%.4f"%(z0+dz0), "%.4f"%(x0+dx00),h_grid,
              "%.4f"%(z1+dz1)]
    L[3,:] = ["%.4f"%(x1+dx11),h_grid,"%.4f"%(z0+dz0), "%.4f"%(x1+dx11),h_grid,
              "%.4f"%(z1+dz1)]
    return L

''' BotTopEdgePartialLines
DESCRIPTION: gets [bottom left, bo right, top le, top re] partial lines at end
             of pixel rows
INPUT: values -> [X[j], X[j]+th_grid, Z[i]+th_grid, Z[i+1], dx0, dx1, dz0, dz1,
                  h_grid]
       * coordinates of 4 bottom points ([x0,z0], [x0,z1], [x1,z1], [x1,z0])
       * dx, dz -> tilt of upper points
       * h_grid -> height of grid
OUTPUT: P -> P_Entries for 4 partial lines at edges of grid,
        parLP -> pointers to [bl, br, tl, tr] partial lines '''
def BotTopEdgePartialLines(entry_index, nextline, values):
    P = []  # output
    pointFlag = True
    ''' additional points: '''
    Points = BotTopRhoPartialPointParams(values[0],values[1],values[3],
                                         values[4],values[5],values[7],
                                         values[8])
    
    for i in range(0,4):  # write Points as P_Entries:
        P = P + P_Entry(entry_index, nextline, 116, Points[i], pointFlag)
        nextline = int(P[-1])
        entry_index = entry_index + 2
    
    parLP = np.zeros(4, dtype=int)  # pointers for lines from entry_index
    parL = BotTopRhoPartialLineParams(values[0],values[1],values[2],values[3],
                                      values[4],values[5],values[6],values[7],
                                      values[8])
    
    for i in range(0,4):  # write parL as P_Entries:
        parLP[i] = entry_index
        P = P + P_Entry(entry_index, nextline, 110, parL[i], pointFlag)
        nextline = int(P[-1])
        entry_index = entry_index + 2
    return [P, parLP]

''' BotTopRhoBoundaries
DESCRIPTION: gets two boundaries of [bottom, top] surfaces between pixel rows
INPUT: rhoE -> pointers to left/right edge lines of surfaces between pixel rows
       rho F/B b/t -> pointers to collection of front/back bottom/top lines of
                      space in between pixel rows
       SP -> pointers to unbounded [bottom, top] surfaces
OUTPUT: 2 boundaries as P_Entries '''
def BotTopRhoBoundaries(rhoE, rhoFb, rhoBb, rhoFt, rhoBt, SP, entry_index,
                        nextline):
    # bottom:
    rhos = [rhoE[0], rhoE[1]]
    P = OuterEdgesRLboundary(rhos, rhoFb, rhoBb, SP[0], entry_index, nextline)
    nextline = int(P[-1])
    entry_index = entry_index + 2
    # top:
    rhos = [rhoE[2], rhoE[3]]
    return P + OuterEdgesRLboundary(rhos, rhoFt, rhoBt, SP[1], entry_index,
                                    nextline)

''' BotTopRhoSurfaces
DESCRIPTION: gets all P_entries of [bottom, top] surfaces between pixels
INPUT: rhoE -> pointers to [left bottom, right bottom, l top, r top] lines of
               surface in between pixel rows
       rho F/B b/t -> front/back bottom/top collection of partial lines along
                      surface in between pixel rows (ordered, left to right)
       btSP -> pointers to unbounded bottom/ top surfaces
OUTPUT: boundaries and bounded surfaces as P_Entries for both surfaces '''
def BotTopRhoSurfaces(rhoE, rhoFb, rhoBb, rhoFt, rhoBt, SP, entry_index,
                      nextline):
    ''' boundaries: '''
    BP = [entry_index, entry_index+2]  # pointers for boundaries
    P = BotTopRhoBoundaries(rhoE, rhoFb, rhoBb, rhoFt, rhoBt, SP,
                            entry_index, nextline)
    nextline = int(P[-1])
    entry_index = entry_index + 4  # + 2*(2 boundaries)
    ''' bounded surfaces '''
    return P + BoundedSurfaces(SP, BP, entry_index, nextline)

'''---------------------------------------------------------------------------
SECTION IV: MAIN ROUTINE '''

''' allP_Entries
creates entire P section as list of entries
OUTPUT:
P, alternating:
    P[0] -> first nextline = 1
    P[1+3i] -> string as written in P section
    P[2+3i] -> type as integer
    P[3+3i] -> line index for NEXT entry ("nextline") '''
def allP_Entries(Angles, GR, th_grid):
    ''' INITIALIZATION: '''
    entry_index = 1  # entry P, stays same value for each element
    # (even if covering several lines)
    # can be used as pointer for this entry
    nextline = 1  # line counter (continually counting used up lines)
    
    alphas = Angles[0]  # per output of determineAllAnglesAndCoordinates
    rhos = Angles[1]
    X = Angles[2]  # positions corresponding to alphas
    Z = Angles[3]  # positions corresponding to rhos
    N_P = Angles[4]  # number of pixels on tile [x,z]
    
    h_grid = th_grid*GR  # height of grid
    
    a = len(alphas)
    r = len(rhos)
    
    P = [1]  # entries list -> output
    # initialize with 1 so that in createIGESfile D section, can iterate linesp
    gapRho = False
    gapAlpha = False
    
    ''' OUTER EDGE LINES (I): '''  # surfaces later (II)
    dx00 = getD(alphas[0], h_grid)  # tilt of left outer edge
    dx11 = getD(alphas[-1], h_grid)  # tilt of right outer edge
    dz00 = getD(rhos[0], h_grid)  # tilt of front outer edge
    dz11 = getD(rhos[-1], h_grid)  # tilt of back outer edge
    
    values = [X[0], X[-1]+th_grid, Z[0], Z[-1]+th_grid, dx00, dx11, dz00, dz11,
              h_grid]
    [P1, oeL] = OuterEdgesLines(entry_index, nextline, values)
    # oeL: pointers to 6 lines of outer edges (right/ left edges excluded)
    P = P + P1
    p = int(len(P1)/3)  # number of additional entries
    nextline = int(P[-1])
    entry_index = entry_index + p*2
    
    ''' BOTTOM & TOP RULED SURFACE: '''
    LBottom = [17,19]  # pointers to 2 antiparallel lines on bottom surface
    LTop = [29,31]  # pointers to 2 antiparallel lines on top surface
    # outer edges are first entries in P section:
        # 8 points + 4 bottom lines + 4 vertical lines + 4 top lines
        # only odd numbers are used
    bSP = entry_index  # pointer to bottom surface
    P = P + RuledSurface(LBottom, entry_index, nextline)
    nextline = int(P[-1])
    entry_index = entry_index + 2
    
    btSP = [bSP, entry_index]  # pointers to [bottom, top] surfaces
    P = P + RuledSurface(LTop, entry_index, nextline)
    nextline = int(P[-1])
    entry_index = entry_index + 2
    
    ''' INNER PIXELS + BOTTOM & TOP BOUNDED SURFACES '''
    dx0 = dx00  # tilt of first alpha leaf = tilt of left outer edge
    dz0 = dz00  # tilt of first rho leaf = tilt of front outer edge
    dxgap = 0  # for alpha gaps, need three dxs simultaneously
    dzgap = 0  # for rho gaps, need three dzs simultaneously
    
    rhoFb = oeL[0]  # pointers to current front edge(s) of rho bottom face
    rhoFt = oeL[6]  # pointers to current front edge(s) of rho top face
    edgeRL = np.zeros(0, int)  # pointers to partial l/r edge lines of grid
    for i in range(0, r-1):  # z
        # start with left edge of pixel, include right edge in loop
        # r-1 -> need to exclude last angle (= right edge) from looping        
        rhoBb = np.zeros(0, int)  # pointers to current back edge(s)
        rhoBt = np.zeros(0, int)  # bottom/ top
        # = pointers to front edges of pixels and lines in between
        nextrhoFb = np.zeros(0, int)  # pointers to next front edge(s)
        nextrhoFt = np.zeros(0, int)  # bottom/ top
        # = pointers to back lines of pixels and lines in between
        if (i+1)%(N_P[1]+1)==0:  # if gap
            dzgap = getD(rhos[i+1], h_grid)
            gapRho = True
            # do nothing
        else:
            dz1 = getD(rhos[i+1], h_grid)  # tilt of right edge
            ''' right/ left edges of rho leaf between pixel rows: '''
            values = [X[0], X[-1]+th_grid, Z[i], Z[i]+th_grid, dx00, dx11, dz0,
                      dz0, h_grid]  # on leaf, dz0 = dz1
            if gapRho:  # if last was gap
                values = [X[0], X[-1]+th_grid, Z[i-1], Z[i]+th_grid, dx00,
                          dx11, dz0, dzgap, h_grid]  # for gap, dz0 ~= dz1
            [P1, rhoE] = BotTopEdgePartialLines(entry_index, nextline, values)
            # rhoE: pointers to [bl, br, tl, tr] lines of surface between pixel
            #       rows
            P = P + P1
            p = int(len(P1)/3)  # number of additional entries
            nextline = int(P[-1])
            entry_index = entry_index + p*2
            edgeRL = np.concatenate((edgeRL, rhoE))
            ''' right/ left lines at end of pixel rows: '''
            values = [X[0], X[-1]+th_grid, Z[i]+th_grid, Z[i+1], dx00, dx11,
                      dz0, dz1, h_grid]
            if gapRho:  # if last was gap
                gapRho = False
                values = [X[0], X[-1]+th_grid, Z[i]+th_grid, Z[i+1], dx00,
                          dx11, dzgap, dz1, h_grid]
            [P1, rhoL] = BotTopEdgePartialLines(entry_index, nextline, values)
            # rhoL: pointers to [bottom left, bo right, top le, top ri] lines
            #       of surface at end of pixel row
            P = P + P1
            p = int(len(P1)/3)  # number of additional entries
            nextline = int(P[-1])
            entry_index = entry_index + p*2
            edgeRL = np.concatenate((edgeRL, rhoL))
            alpLb = rhoL[0]
            alpLt = rhoL[2]
            
            for j in range(0, a-1):  # x
                if (j+1)%(N_P[0]+1)==0:  # if gap
                    dxgap = dx0
                    dx0 = getD(alphas[j+1], h_grid)
                    gapAlpha = True
                    # do nothing
                else:
                    dx1 = getD(alphas[j+1], h_grid)
                    ''' pixel boardering surfaces: '''
                    values = [X[j]+th_grid, X[j+1], Z[i]+th_grid, Z[i+1], dx0,
                              dx1, dz0, dz1, h_grid]
                    [P1, L1] = P_EntriesOnePixel(entry_index, nextline, values)
                    P = P + P1
                    p = int(len(P1)/3)  # number of additional entries
                    nextline = int(P[-1])
                    entry_index = entry_index + p*2
                    ''' (alpha) bottom and top faces between pixels: '''
                    alpRb = L1[3]
                    alpRt = L1[11]
                    
                    values = [X[j], X[j]+th_grid, Z[i]+th_grid, Z[i+1], dx0,
                              dx0, dz0, dz1, h_grid]  # dx1 = dx0 if on leaf 
                    if gapAlpha:  # if last was gap
                        gapAlpha = False
                        values = [X[j-1], X[j]+th_grid, Z[i]+th_grid, Z[i+1],
                                  dxgap, dx0, dz0, dz1, h_grid]  # dx1 ~= dx0
                    alps = [alpLb, alpRb, alpLt, alpRt]
                    # pointers to lines in pixels that bind edge faces
                    
                    [P1, pL] = BotTopAlphaSurfaces(alps, btSP, entry_index,
                                                   nextline, values)
                    P = P + P1
                    p = int(len(P1)/3)  # number of additional entries
                    nextline = int(P[-1])
                    entry_index = entry_index + p*2
                    
                    rhoBb = np.concatenate((rhoBb, [pL[0],L1[0]]))
                    rhoBt = np.concatenate((rhoBt, [pL[2],L1[8]]))
                    nextrhoFb = np.concatenate((nextrhoFb, [pL[1], L1[2]]))
                    nextrhoFt = np.concatenate((nextrhoFt, [pL[3], L1[10]]))
                    
                    alpLb = L1[1]  # left bottom edge of next face
                    alpLt = L1[9]  # left top edge of next face
                    dx0 = dx1
                    # tilt of left edge, reuse for right edge of next pixel
            ''' last alpha leaf bottom and top face: '''
            alpRb = rhoL[1]
            alpRt = rhoL[3]
            alps = [alpLb, alpRb, alpLt, alpRt]
            
            values = [X[j+1], X[j+1]+th_grid, Z[i]+th_grid, Z[i+1], dx0, dx0,
                              dz0, dz1, h_grid]
            [P1, pL] = BotTopAlphaSurfaces(alps, btSP, entry_index, nextline,
                                           values)
            P = P + P1
            p = int(len(P1)/3)  # number of additional entries
            nextline = int(P[-1])
            entry_index = entry_index + p*2
            
            rhoBb = np.concatenate((rhoBb, [pL[0]]))
            rhoBt = np.concatenate((rhoBt, [pL[2]]))
            nextrhoFb = np.concatenate((nextrhoFb, [pL[1]]))
            nextrhoFt = np.concatenate((nextrhoFt, [pL[3]]))
            ''' rho leaf bottom and top face between pixel rows: '''
            P1 = BotTopRhoSurfaces(rhoE, rhoFb, rhoBb, rhoFt, rhoBt, btSP,
                                   entry_index, nextline)
            P = P + P1
            p = int(len(P1)/3)  # number of additional entries
            nextline = int(P[-1])
            entry_index = entry_index + p*2
            
            rhoFb = nextrhoFb
            rhoFt = nextrhoFt
            
            dz0 = dz1  # tilt of back edge, reuse for front edge of next pixel
            dx0 = dx00  # reset to first leaf
    ''' last rho leaf bottom and top face: '''
    values = [X[0], X[-1]+th_grid, Z[-1], Z[-1]+th_grid, dx00, dx11, dz0, dz0,
              h_grid]
    [P1, rhoE] = BotTopEdgePartialLines(entry_index, nextline, values)
    # rhoE: pointers to [bl, br, tl, tr] lines of surface between pixel
    #       rows
    P = P + P1
    p = int(len(P1)/3)  # number of additional entries
    nextline = int(P[-1])
    entry_index = entry_index + p*2
    edgeRL = np.concatenate((edgeRL, rhoE))
    
    rhoBb = oeL[1]  # pointer to bottom back edge of grid
    rhoBt = oeL[7]  # pointer to top back edge of grid
    P1 = BotTopRhoSurfaces(rhoE, rhoFb, rhoBb, rhoFt, rhoBt, btSP, entry_index,
                           nextline)
    P = P + P1
    p = int(len(P1)/3)  # number of additional entries
    nextline = int(P[-1])
    entry_index = entry_index + p*2
    ''' OUTER EDGE SURFACES (II): '''
    return P + OuterEdgesSurfaces(oeL, edgeRL, btSP, entry_index, nextline)