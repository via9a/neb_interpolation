import sys
import os
import numpy as np

"""
Author: Vilhjalmur Asgeirsson (UI, 2018)
"""


def Displacement(ndim, nim, R):

    displ = np.zeros(shape=(nim,))
    for i in range(1, nim):
        R0 = R[(i - 1) * ndim:(i) * ndim]
        R1 = R[(i) * ndim:(i + 1) * ndim]
        dR = R1-R0
        displ[i] = displ[i - 1] + np.sqrt(np.dot(dR.T,dR))

    return displ

def LinearInterpolateData(nlen, xData, yData, xnew):

    i = 0
    if (xnew >= xData[nlen - 2]):
        i = nlen - 2
    else:
        while (xnew > xData[i + 1]):
            i = i + 1

    xL = xData[i]
    yL = yData[i]
    xR = xData[i + 1]
    yR = yData[i + 1]

    if (xnew < xL):
        yR = yL
    if (xnew > xR):
        yL = yR
    dydx = (yR - yL) / (xR - xL)

    return yL + dydx * (xnew - xL)


def GenerateNewPath(ndim, nim, npoints, S, R):

    newR = np.zeros(shape=(ndim * npoints, 1))
    xi = np.linspace(S[0], S[-1], npoints)

    for i in range(ndim):
        Rdof = np.zeros(shape=(nim, 1))
        dRdof = np.zeros(shape=(nim, 1))
        for j in range(nim):
            Rdof[j] = R[j * ndim + i]

        for xpt in range(npoints):
            new_y = LinearInterpolateData(nim, S, Rdof, xi[xpt])
            newR[xpt * ndim + i] = new_y

    return newR

def ReadFirstLineOfFile(fname):
    with open(fname) as f:
        first_line = f.readline()
    return first_line

def WriteTraj(fname, ndim, nim, R, E, symb3):

    if len(symb3) != ndim:
        raise RuntimeError("Error in WriteTraj. Dimension mismatch: symb3")
    if len(E) != nim:
        raise RuntimeError("Error in WriteTraj. Dimension mismatch: E")
    if len(R) != ndim*nim:
        raise RuntimeError("Error in WriteTraj. Dimension mismatch: R")
    
    f = open(fname, 'w')
    natoms = ndim / 3
    for i in range(nim):
        hnit = R[i * ndim:(i + 1) * ndim]
        f.write("%i \n" % natoms )
        f.write('E=%12.8f \n' % E[i])
        for j in range(0, len(hnit), 3):
            f.write('%s %12.8f %12.8f %12.8f\n' % (symb3[j].strip(), hnit[j], hnit[j + 1], hnit[j + 2]))
    f.close()
    
    return None

def ReadTraj(fname):
    if not os.path.isfile(fname):
        raise RuntimeError("File %s not found " % fname)

    extension = fname.split('.')[-1]
    
    # Get first line of the file (i.e. number of atoms)
    first_line = ReadFirstLineOfFile(fname)

    try:
        natoms = int(first_line)
    except:
        raise TypeError("%s is not correctly formatted trajectory file")
    
    ndim = natoms * 3

    # begin by reading the contents of the file to a list
    contents = []
    f = open(fname).readlines()
    for i, line in enumerate(f):
        contents.append(line)

    # get number of lines and hence number of images
    number_of_lines = i + 1
    nim = int((number_of_lines) / (natoms + 2))

    RPATH = []
    ind = 0
    for i in range(nim):
        symb3 = []
        ind = ind + 2
        for j in range(natoms):
            geom_line = contents[ind]
            geom_line = geom_line.split()
            symb3.append(geom_line[0].strip())
            symb3.append(geom_line[0].strip())
            symb3.append(geom_line[0].strip())
            RPATH.append(float(geom_line[1]))
            RPATH.append(float(geom_line[2]))
            RPATH.append(float(geom_line[3]))
            ind += 1

    RPATH = np.reshape(RPATH, (nim * ndim, 1))

    return RPATH, ndim, nim, symb3


if __name__ == "__main__":

    """
    Script to extend a trajectory from ORCA NEB run to include more/less images
    by using simple linear interpolation. New number of images is given by 
    the variable npts

    requires: numpy

    Usage: python neb_interpolate_path.py basename_MEP.trj npts<int>
    (in the given order)

    Authors: Vilhjalmur Asgeirsson (UI, 2018)
    """

    # ============================================
    # default values
    # ============================================
    
    fname = 'orca_MEP.trj'
    npts = 100

    # ============================================
    # Get inp arguments
    # ============================================
    
    for i in range(len(sys.argv)):
        if i == 1:
            fname = sys.argv[1]
        if i == 2:
            try:
                npts  = int(sys.argv[2])
            except:
                raise TypeError("Int. expected as a second arg")
        
    # ============================================
    # Read trajectory file
    # ============================================
    basename = fname.split('.')[0]
    fname_output = basename+'_extended.xyz'
    R, ndim, nim, symb3 = ReadTraj(fname)
    
    # ============================================
    # Perform interpolation
    # ============================================
    S = Displacement(ndim, nim, R)
    newR = GenerateNewPath(ndim, nim, npts, S, R)
    
    # ============================================
    # Write new trajectory file
    # ============================================
    E = np.zeros(shape=(npts,1))
    WriteTraj(fname_output, ndim, npts, newR, E, symb3)

    # Done.
