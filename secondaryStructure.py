import prody as dy
import matplotlib.pyplot as plt
import os

import evaluateinterface


path = "/home/glenn/Documents/WebNma/best"
path = "/home/glenn/work/benchmark5_attract/"
names = []
#file =open(os.path.join(path, "dir.txt"))

# for line in file.readlines():
#     names.append(line.strip("\n"))
#     print line
# file.close()

names =['1GCQ']

names =[]

with open("/home/glenn/work/benchmark5_attract/protlist") as f:
    for line in f.readlines():
        names.append(line.split()[0])


for name in names:
    protPath = os.path.join(path, name)
    proteinName = ["{}A-unbound".format(name),"{}B-unbound".format(name) ]
    for prot in proteinName:
        inPDB = os.path.join(protPath, prot+".pdb")
        print inPDB
        outSTRIDE = os.path.join(protPath, prot)
        outSecondary = os.path.join(protPath,prot+".sec")
        dy.execSTRIDE(inPDB, outSTRIDE)
        pdb = dy.parsePDB(inPDB)
        dy.parseSTRIDE(outSTRIDE+".stride", pdb)
        secondaryStruc = pdb.getSecstrs()
        resIdx = pdb.getResindices()
        secFile = open(outSecondary, 'w+')
        curres = -1
        secFile.write("residueNumber\n".format(sec))
        for sec,res in zip(secondaryStruc, resIdx):
            if res != curres:
                secFile.write("{}\n".format(sec))
                curres = res
        secFile.close()



def createSecStruc(filename, outSTRIDE):
    dy.execSTRIDE(filename, outSTRIDE)
    pdb = dy.parsePDB(filename)
    dy.parseSTRIDE(outSTRIDE + ".stride", pdb)
    secondaryStruc = pdb.getSecstrs()
   # resIdx = pdb.getResindices()
    return secondaryStruc


def read_modes(filename_modes):

    with open(filename_modes, 'r') as f:
        lines = f.readlines()
        d_modes = {}
        mode_idx = 1
        for i, line in enumerate(lines):
            if len(line.split()) == 2 and int(line.split()[0]) == mode_idx:
                x = []
                y = []
                z = []
                d_modes[mode_idx ] = ( float(line.split()[1]), x,y,z )
                mode_idx += 1
                #eigenvalues[mode_idx - 2] = float(line.split()[1])
            elif len(line.split()) == 3:
                d_modes[mode_idx - 1][1].append( float(line.split()[0]) )
                d_modes[mode_idx - 1][2].append( float(line.split()[1]) )
                d_modes[mode_idx - 1][3].append( float(line.split()[2]) )

    return d_modes

def fluctuation(x,y,z):
    return x*x+ y*y+z*z

def reduceToResidueMode(mx, my, mz, length):
    currline = [10000,10000,10000]
    mrx, mry, mrz = [],[],[]
    for i in range(length):
        if(mx[i] != currline[0] or my[i] != currline[1] or mz[i] != currline[2] ):
            currline = [mx[i] ,my[i] ,mz[i] ]
            mrx.append(mx[i])
            mry.append(my[i])
            mrz.append(mz[i])

    return np.asarray(mrx),np.asarray(mry),np.asarray(mrz)

import pandas as pd

def plotFulctuation( path , protein, numModes, hinsen = False):
    protPath = os.path.join(path, protein)
    fig, axes = plt.subplots(2, 2)
    if (hinsen):
        fig, axes = plt.subplots(3,2)
    inSec = os.path.join(protPath, protein +"A-unbound.sec")
    filename_modes = os.path.join(protPath, 'input/'+protein + "-receptor-for-docking-10-modes.dat")

    df = pd.read_csv(inSec)

    axes[0,0].plot(df['residueNumber'])
    dm = read_modes(filename_modes)
    ylim = 0
    for i in range(numModes):
        x, y, z = reduceToResidueMode(dm[i+1][1], dm[i+1][2], dm[i+1][3], len(dm[i+1][1]))
        fluc = fluctuation(x, y, z)
        fluc = np.sqrt(fluc)
        axes[1,0].plot(fluc, label = "mode {}".format(i+1))
       # axes[2, 0].plot(np.fft.fft(fluc), label="mode {}".format(i + 1))
        ylim = max(ylim, max(fluc[:len(fluc) - 1]))
        axes[1,0].set_ylim(0,ylim)

    if (hinsen):
        filename_modes = os.path.join(protPath, 'input/' + protein + "-receptor-for-docking-10-hinsen.dat")

        dm = read_modes(filename_modes)
        ylim = 0
        for i in range(numModes):
            x, y, z = reduceToResidueMode(dm[i + 1][1], dm[i + 1][2], dm[i + 1][3], len(dm[i + 1][1]))
            fluc = fluctuation(x, y, z)
            fluc = np.sqrt(fluc)
            axes[2, 0].plot(fluc, label="hinsen {}".format(i + 1))
            ylim = max(ylim, max(fluc[:len(fluc) - 1]))
            axes[2, 0].set_ylim(0, ylim)

    inSec = os.path.join(protPath, protein +"B-unbound.sec")
    filename_modes = os.path.join(protPath, 'input/'+protein + "-ligand-for-docking-10-modes.dat")

    df = pd.read_csv(inSec)

    axes[0,1].plot(df['residueNumber'])
    dm = read_modes(filename_modes)
    ylim = 0
    for i in range(numModes):
        x, y, z = reduceToResidueMode(dm[i+1][1], dm[i+1][2], dm[i+1][3], len(dm[i+1][1]))
        fluc = fluctuation(x, y, z)
        fluc = np.sqrt(fluc)
        axes[1,1].plot(fluc, label = "mode {}".format(i+1))
       # axes[2, 1].plot(np.fft.fft(fluc)**2, label="mode {}".format(i + 1))
        ylim = max(ylim,max(fluc[:len(fluc)-1]))
        axes[1,1].set_ylim(0,ylim)
    axes[1, 1].legend()
    axes[1, 0].legend()

    if (hinsen):
        filename_modes = os.path.join(protPath, 'input/'+protein + "-ligand-for-docking-10-hinsen.dat")

        dm = read_modes(filename_modes)
        ylim = 0
        for i in range(numModes):
            x, y, z = reduceToResidueMode(dm[i+1][1], dm[i+1][2], dm[i+1][3], len(dm[i+1][1]))
            fluc = fluctuation(x, y, z)
            fluc = np.sqrt(fluc)
            axes[2,1].plot(fluc, label = "hinsen {}".format(i+1))
            ylim = max(ylim,max(fluc[:len(fluc)-1]))
            axes[2,1].set_ylim(0,ylim)
    axes[1, 1].legend()
    axes[1, 0].legend()
    if (hinsen):
        axes[2, 1].legend()
        axes[2, 0].legend()

plotFulctuation( path , '1AVX', 5)


names=['3MXW']
for name in names:
    proteinName = ["{}A-unbound".format(name), "{}B-unbound".format(name)]
    for prot in proteinName:
        fileSecStruc = "{}/{}/{}.sec".format(path,name,prot)
        fileDeformation = "{}/{}/webnma/{}_fluctuationsplot.dat".format(path,name,prot)
        secFile = open(fileSecStruc)
        sec =[]
        for line in secFile.readlines():
            sec.append(line.split()[0])
        secFile.close()

        defFile = open(fileDeformation)
        defo = []
        for line in defFile.readlines():
            defo.append(float(line.split()[1]))
        plt.clf()
        plt.plot(sec)
        plt.plot(defo)
        plt.show()