import os
#from  load_pdbs import ProteinEsemble
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



class Path:
    def __init__(self, path):
        self.path = path
    def filename(self, filename):
        return os.path.join( self.path, filename )
    def isfile(self, filename):
        return os.path.isfile( os.path.join( self.path, filename))


def isFile(filename):
    return os.path.isfile(filename)

def readLineFile(filename, colIdx, sep = " ", type = float):
    try:
        with open(filename) as f:
            lines = f.readlines()
            data = np.zeros(len(lines))
            #data=list[len(lines)]
            for i, line in enumerate(lines):
                #data[i] = np.fromstring(line, dtype=type, sep=sep)[colIdx]
                data[i] = type(line.split()[colIdx])
        return data
    except:
        return -1

def readRMSD(filename):
    return readLineFile(filename, 1)

def readFnat(filename):
    return readLineFile(filename, 0)



def getEnergyfromFile( filename):
    energies = []
    try:
        with open(filename, 'r') as f:
            count = 0
            idx = 1
            for line in f.readlines():
                if line.startswith("#{}".format(idx)):
                    number =  float(line[1])
                    idx += 1
                    count += 1
                elif line.startswith("## Energy:"):
                   # print line
                    energy = float(line.split()[2])
                    count += 1
                if count == 2:
                    energies.append( energy)
                    count = 0
        if idx == 1:
            energies = [-100000]
        return np.asarray( energies )
    except:
        print("execpt get energy")
        return np.asarray( [-100000] )

#testcase eread
#fn = '/home/glenn/work/benchmark5_best/1IJK/dG_mr5_ml1_ev0p1_sO_c50_mr5_ml1_ev0p1_hinsen/1IJK-receptor-for-docking-sorted-dr.dat'
#energies = getEnergyfromFile( fn)

def getCapriRanking(Irmsd, fnat):
        if fnat >= 0.5 and Irmsd <= 1.0:
            return 3
        elif fnat >= 0.3 and Irmsd <= 2.0:
            return 2
        elif fnat >= 0.1 and fnat < 0.3 and Irmsd < 4:
            return 1
        else:
            return 0


def getCapriList(IrmsdData, fnatData):
    list_capri = np.zeros(len(IrmsdData))
    for i, (irmsd, fnat) in enumerate(zip(IrmsdData, fnatData)):
        list_capri[i] = getCapriRanking(irmsd, fnat)
    return list_capri


def capriCount(list_capri):
    count_1 = 0
    count_2 = 0
    count_3 = 0
    for rank in list_capri:
        if rank == 1:
            count_1 += 1
        elif rank == 2:
            count_2 += 1
        elif rank == 3:
            count_3 += 1
    return count_1, count_2, count_3

#irmsd = readRMSD("/home/glenn/work/benchmark5_best/1IJK/dG_mr5_ml1_ev0p1_sO_c50_mr5_ml1_ev0p1_hinsen/1IJK-receptor-for-docking-irmsd.dat")
#fnat = readFnat("/home/glenn/work/benchmark5_best/1IJK/dG_mr5_ml1_ev0p1_sO_c50_mr5_ml1_ev0p1_hinsen/1IJK-receptor-for-docking-fnat.dat")

#capriList = getCapriList(irmsd, fnat)
#one_star, two_star, three_star = capriCount(capriList)
#print one_star, two_star, three_star


def getModeEVScale(benchmarkname):
    index_start = benchmarkname.find("ev") + 2
    index_dec = benchmarkname[index_start:].find("p") + 1 + index_start
    endIdx = benchmarkname[index_dec:].find('_') + index_dec

    return float(benchmarkname[index_start:index_dec - 1]) + float("0." + benchmarkname[index_dec:endIdx])

#print getModeEVScale(bm)

def getNumModesFromString(benchmarkname, protName):
    modeMarker = "m"+protName
    index_modesLig = benchmarkname.find(modeMarker) + 2
    index_end = benchmarkname[index_modesLig:].find('_') + index_modesLig
    return int(benchmarkname[index_modesLig:index_end])

#print getNumModesFromString(bm, 'l')

#load benchmark from list of names and benchmarks and path
#classify ( number of modes and scale of the eigenvalues
#

def getPath(basePath,benchmarkName, proteinName):
    return os.path.join(basePath, os.path.join(proteinName, benchmarkName))


fext = {'rmsd':"-rmsd.result",
     'irmsd': "-irmsd.dat",
     'fnat' : "-fnat.dat",
     'energy' : "-sorted-dr.dat",
     'pdb' : "-receptor-for-docking"}



def readStringList(filename, colIdx):
    with open(filename) as f:
        lines = f.readlines()
        list = []
        for line in lines:
            list.append(line.split()[colIdx])
        return  list

def evalData(protList, benchmarkList,basePath, fileout):
    dataList=[]
    for bm in benchmarkList:
        print bm
        for protein in protList:


            bmPath = Path(getPath(basePath=basePath, benchmarkName=bm, proteinName= protein))
            evScale = "{:2.6f}".format(getModeEVScale(bm))
            numModesRec = getNumModesFromString(bm, 'r')
            numModesLig = getNumModesFromString(bm, 'l')
            baseFile = protein+fext['pdb']

            # energyFile = baseFile + fext['energy']
            # energies = []
            # if(bmPath.isfile(energyFile)):
            #     energies = getEnergyfromFile(bmPath.filename(energyFile))
            # rmsdFile = baseFile + fext['rmsd']
            # rmsd =[]
            # if (bmPath.isfile(rmsdFile)):
            #     rmsd = np.asarray(readRMSD(bmPath.filename(rmsdFile)))
            irmsdFile = baseFile + fext['irmsd']
            irmsd = []
            if (bmPath.isfile(irmsdFile)):
                irmsd = np.asarray(readRMSD(bmPath.filename(irmsdFile)))
            fnatFile = baseFile + fext['fnat']
            fnat =[]
            if (bmPath.isfile(fnatFile)):
                fnat = np.asarray(readFnat(bmPath.filename(fnatFile)))
            #print fnat
            if type(fnat) == int or type(irmsd) == int  :
                print "file not available for ", protein
                continue
            if len(fnat) > 1 and len(irmsd) > 1:
                capriList = getCapriList(irmsd, fnat)
                one_star, two_star, three_star = capriCount(capriList[:50])
                one_star10, two_star10, three_star10 = capriCount(capriList[:10])
                one_10, two_10, three_10= capriCount(capriList[:10])
                one_50, two_50, three_50 = capriCount(capriList[:50])
                if(len(capriList) >50):
                    one_100, two_100, three_100 = capriCount(capriList[:100])
                    one_1000, two_1000, three_1000 = capriCount(capriList[:1000])
                else:
                    one_100=0
                    two_100=0
                    three_100=0
                    one_1000 = 0
                    two_1000 = 0
                    three_1000 = 0
                #print bm.find('hinsen')
                if bm.find('hinsen') != -1:
                    brownian = 1
                else:
                    brownian = 0
                if bm.find('omodes') != -1:
                    omodes = 1
                else:
                    omodes = 0
                if bm.find('bound') != -1:
                    bound = 1
                else:
                    bound = 0
                row = {'protein':protein,'numModesRec':numModesRec,'numModesLig':numModesLig,'evScale':evScale, 'one_star_50':one_star,'two_star_50':two_star,'three_star_50':three_star,'one_star_10':one_star10,'two_star_10':two_star10,'three_star_10':three_star10, 'brownian': brownian,'omodes':omodes, 'bound':bound,'one_10':one_10>0, 'two_10':two_10>0, 'three_10':three_10>0
                    , 'one_50': one_50>0, 'two_50': two_50>0, 'three_50': three_50>0,'one_100':one_100>0, 'two_100':two_100>0, 'three_100':three_100>0,'one_1000':one_1000>0, 'two_1000':two_1000>0, 'three_1000':three_1000>0}
                dataList.append(row)
                #print protein,numModesLig, numModesRec, evScale,one_star, two_star, three_star, brownian
            else:
                print "false", protein
    frame = pd.DataFrame(dataList)
    frame.to_csv(fileout)

def saveRmsdEnergy(protList, benchmarkList,basePath,outbasepath):
    for protein in protList:
        for bm in benchmarkList:
            bmPath = Path(getPath(basePath=basePath, benchmarkName=bm, proteinName=protein))
            print getPath(basePath=basePath, benchmarkName=bm, proteinName=protein)
            evScale = getModeEVScale(bm)
            numModesRec = getNumModesFromString(bm, 'r')
            numModesLig = getNumModesFromString(bm, 'l')
            baseFile = protein + fext['pdb']

            energyFile = baseFile + fext['energy']
            energies = []
            if (bmPath.isfile(energyFile)):
                energies = getEnergyfromFile(bmPath.filename(energyFile))
            rmsdFile = baseFile + fext['rmsd']
            rmsd = []
            if (bmPath.isfile(rmsdFile)):
                rmsd = np.asarray(readRMSD(bmPath.filename(rmsdFile)))
            if  type(rmsd) == int:
                print "file not available for ", protein
                continue
            option = ""
            if bm.find('hinsen') != -1:
                option += "-b"
            if bm.find('omodes') != -1:
                option += "-om"
            outPath = "{}/prot{}-mr{}-ml{}-ev{}-{}.rmsdplot".format(outbasepath,protein,numModesRec, numModesLig,evScale,option)


            if len(rmsd) >1 and len(energies)> 1:
                frame = pd.DataFrame()
                frame['rmsd'] = rmsd
                frame['energy'] = energies
                frame.to_csv(outPath)

def plotRmsdEnergy(protList, benchmarkList,basePath):
    for protein in protList:
        for bm in benchmarkList:
            bmPath = Path(getPath(basePath=basePath, benchmarkName=bm, proteinName=protein))
            print getPath(basePath=basePath, benchmarkName=bm, proteinName=protein)
            evScale = getModeEVScale(bm)
            numModesRec = getNumModesFromString(bm, 'r')
            numModesLig = getNumModesFromString(bm, 'l')
            baseFile = protein + fext['pdb']

            energyFile = baseFile + fext['energy']
            energies = []
            if (bmPath.isfile(energyFile)):
                energies = getEnergyfromFile(bmPath.filename(energyFile))
            rmsdFile = baseFile + fext['rmsd']
            rmsd = []
            if (bmPath.isfile(rmsdFile)):
                rmsd = np.asarray(readRMSD(bmPath.filename(rmsdFile)))


            if len(rmsd) >1 and len(energies)> 1:
                plt.scatter(rmsd, energies, alpha= 0.4, label = bm)
                plt.legend()
                #xy = np.vstack([rmsd, energies])
                #z = gaussian_kde(xy)(xy)
                #fig, ax = plt.subplots()
                #ax.scatter(rmsd, energies, c=z, s=100, edgecolor='', alpha = 0.5)


def getEnergyDelta(protList, benchmarkList,basePath,numtop):
    d = {}
    for protein in protList:
        sum = []
        yep = True
        for bm in benchmarkList:
            bmPath = Path(getPath(basePath=basePath, benchmarkName=bm, proteinName=protein))
            #print getPath(basePath=basePath, benchmarkName=bm, proteinName=protein)
            evScale = getModeEVScale(bm)
            numModesRec = getNumModesFromString(bm, 'r')
            numModesLig = getNumModesFromString(bm, 'l')
            baseFile = protein + fext['pdb']

            energyFile = baseFile + fext['energy']

            if (bmPath.isfile(energyFile)):
                energies = []
                energies = getEnergyfromFile(bmPath.filename(energyFile))

                if len(energies)>= numtop:
                    energies = energies[:numtop]
                    sum.append(np.sum(energies))
                else:
                    yep = False
            else:
                yep = False


        if yep is True:
            delta = sum[0]- sum[1]
            delta /= numtop
            d[protein] = delta
            print(protein , delta)
    return d


def loadRmsd(protList, benchmarkList,basePath):
    dict={}
    for protein in protList:
        for bm in benchmarkList:
            bmPath = Path(getPath(basePath=basePath, benchmarkName=bm, proteinName=protein))

            baseFile = protein + fext['pdb']
            rmsdFile = baseFile + fext['rmsd']
            rmsd = []
            if (bmPath.isfile(rmsdFile)):
                rmsd  = readRMSD(bmPath.filename(rmsdFile))
                if len(rmsd) > 1:
                    dict[protein] = np.asarray(rmsd)
    return dict

def findMinRMSD(protList, benchmarkList,basePath, referencermsd, num, minnum):
    minsum = np.zeros(minnum)
    for protein in protList:
        for bm in benchmarkList:
            bmPath = Path(getPath(basePath=basePath, benchmarkName=bm, proteinName=protein))
            #print getPath(basePath=basePath, benchmarkName=bm, proteinName=protein)
            evScale = getModeEVScale(bm)
            numModesRec = getNumModesFromString(bm, 'r')
            numModesLig = getNumModesFromString(bm, 'l')
            baseFile = protein + fext['pdb']
            rmsdFile = baseFile + fext['rmsd']
            rmsd = []
            if (bmPath.isfile(rmsdFile)):
                rmsd = np.asarray(readRMSD(bmPath.filename(rmsdFile)))[:num]
                if len(rmsd) > 1:
                    sum = np.sort(referencermsd[protein][:num]-rmsd  )[:minnum]
                    if np.sum(sum < minsum)> minnum/2:
                        minModesRec = numModesRec
                        minModesLig = numModesLig
                        minscale = evScale
                        minsum = sum
                        print minsum, numModesRec, numModesLig, protein


def findMinRMSD(protList, benchmarkList,basePath, referencermsd, num, minnum):
    minsum = np.zeros(minnum)
    for protein in protList:
        for bm in benchmarkList:
            bmPath = Path(getPath(basePath=basePath, benchmarkName=bm, proteinName=protein))
            #print getPath(basePath=basePath, benchmarkName=bm, proteinName=protein)
            evScale = getModeEVScale(bm)
            numModesRec = getNumModesFromString(bm, 'r')
            numModesLig = getNumModesFromString(bm, 'l')
            baseFile = protein + fext['pdb']
            rmsdFile = baseFile + fext['rmsd']
            rmsd = []
            if (bmPath.isfile(rmsdFile)):
                rmsd = np.asarray(readRMSD(bmPath.filename(rmsdFile)))[:num]
                if len(rmsd) > 1:
                    sum = np.sort(referencermsd[protein][:num]-rmsd  )[:minnum]
                    if np.sum(sum < minsum)> minnum/2:
                        minModesRec = numModesRec
                        minModesLig = numModesLig
                        minscale = evScale
                        minsum = sum
                        print minsum, numModesRec, numModesLig, protein

#findMinRMSD(protList, ['dG_mr1_ml1_ev2p0_sO_c50_mr1_ml1_ev2p0_omodes'],"/home/glenn/work/benchmark5_attract", rmsdref, 20, 3)



path = "/home/glenn/work/test_hinsen/"
filename  = "/home/glenn/Documents/Masterarbeit/analysis/largerun/benchmark5_3pc8_hinsen"
protList = readStringList(path + "protlist",0                          )
#benchmarkList = readStringList("/home/glenn/work/benchmark5_best/benchlist",0                          )
#saveRmsdEnergy(protList, benchmarkList,"/home/glenn/work/benchmark5_attract","/home/glenn/Documents/Masterarbeit/analysis/largerun/rmsdplot")
print "eval benchmark"


benchmarkList=[

'dG_mr10_ml10_ev0p00030_sO_c50_mr10_ml10_ev0p00030hinsen',
'dG_mr10_ml10_ev0p00035_sO_c50_mr10_ml10_ev0p00035hinsen',
'dG_mr10_ml10_ev0p00050_sO_c50_mr10_ml10_ev0p00050hinsen',
'dG_mr10_ml10_ev0p00080_sO_c50_mr10_ml10_ev0p00080hinsen',
'dG_mr10_ml10_ev0p70000_sO_c50_mr10_ml10_ev0p70000',
'dG_mr10_ml10_ev1p00000_sO_c50_mr10_ml10_ev1p00000',
'dG_mr10_ml10_ev1p30000_sO_c50_mr10_ml10_ev1p30000',
'dG_mr10_ml10_ev1p70000_sO_c50_mr10_ml10_ev1p70000',
'dG_mr0_ml10_ev0p00005_sO_c50_mr0_ml10_ev0p00005hinsen'  ,
'dG_mr10_ml10_ev2p00000_sO_c50_mr10_ml10_ev2p00000',
'dG_mr0_ml10_ev0p00008_sO_c50_mr0_ml10_ev0p00008hinsen'  ,
'dG_mr10_ml10_ev2p50000_sO_c50_mr10_ml10_ev2p50000',
'dG_mr0_ml10_ev0p00010_sO_c50_mr0_ml10_ev0p00010hinsen'  ,
'dG_mr1_ml0_ev0p00005_sO_c50_mr1_ml0_ev0p00005',
'dG_mr0_ml10_ev0p00012_sO_c50_mr0_ml10_ev0p00012hinsen'  ,
'dG_mr1_ml0_ev0p00005_sO_c50_mr1_ml0_ev0p00005_hinsen',
'dG_mr0_ml10_ev0p00014_sO_c50_mr0_ml10_ev0p00014hinsen'  ,
'dG_mr1_ml0_ev0p00008_sO_c50_mr1_ml0_ev0p00008',
'dG_mr0_ml10_ev0p00015_sO_c50_mr0_ml10_ev0p00015hinsen'  ,
'dG_mr1_ml0_ev0p00008_sO_c50_mr1_ml0_ev0p00008_hinsen',
'dG_mr0_ml10_ev0p00016_sO_c50_mr0_ml10_ev0p00016hinsen'  ,
'dG_mr1_ml0_ev0p00010_sO_c50_mr1_ml0_ev0p00010',
'dG_mr0_ml10_ev0p00017_sO_c50_mr0_ml10_ev0p00017hinsen'  ,
'dG_mr1_ml0_ev0p00010_sO_c50_mr1_ml0_ev0p00010_hinsen',
'dG_mr0_ml10_ev0p00018_sO_c50_mr0_ml10_ev0p00018hinsen'  ,
'dG_mr1_ml0_ev0p00012_sO_c50_mr1_ml0_ev0p00012',
'dG_mr0_ml10_ev0p00020_sO_c50_mr0_ml10_ev0p00020hinsen'  ,
'dG_mr1_ml0_ev0p00012_sO_c50_mr1_ml0_ev0p00012_hinsen',
'dG_mr0_ml10_ev0p00030_sO_c50_mr0_ml10_ev0p00030hinsen'  ,
'dG_mr1_ml0_ev0p00014_sO_c50_mr1_ml0_ev0p00014',
'dG_mr0_ml10_ev0p00035_sO_c50_mr0_ml10_ev0p00035hinsen'  ,
'dG_mr1_ml0_ev0p00014_sO_c50_mr1_ml0_ev0p00014_hinsen',
'dG_mr0_ml10_ev0p00050_sO_c50_mr0_ml10_ev0p00050hinsen'  ,
'dG_mr1_ml0_ev0p00015_sO_c50_mr1_ml0_ev0p00015',
'dG_mr0_ml10_ev0p00080_sO_c50_mr0_ml10_ev0p00080hinsen'  ,
'dG_mr1_ml0_ev0p00015_sO_c50_mr1_ml0_ev0p00015_hinsen',
'dG_mr0_ml10_ev0p70000_sO_c50_mr0_ml10_ev0p70000'  ,
'dG_mr1_ml0_ev0p00016_sO_c50_mr1_ml0_ev0p00016',
'dG_mr0_ml10_ev1p00000_sO_c50_mr0_ml10_ev1p00000'  ,
'dG_mr1_ml0_ev0p00016_sO_c50_mr1_ml0_ev0p00016_hinsen',
'dG_mr0_ml10_ev1p30000_sO_c50_mr0_ml10_ev1p30000'  ,
'dG_mr1_ml0_ev0p00017_sO_c50_mr1_ml0_ev0p00017',
'dG_mr0_ml10_ev1p70000_sO_c50_mr0_ml10_ev1p70000'  ,
'dG_mr1_ml0_ev0p00017_sO_c50_mr1_ml0_ev0p00017_hinsen',
'dG_mr0_ml10_ev2p00000_sO_c50_mr0_ml10_ev2p00000'  ,
'dG_mr1_ml0_ev0p00018_sO_c50_mr1_ml0_ev0p00018',
'dG_mr0_ml10_ev2p50000_sO_c50_mr0_ml10_ev2p50000'  ,
'dG_mr1_ml0_ev0p00018_sO_c50_mr1_ml0_ev0p00018_hinsen',
'dG_mr0_ml1_ev0p00005_sO_c50_mr0_ml1_ev0p00005'  ,
'dG_mr1_ml0_ev0p00020_sO_c50_mr1_ml0_ev0p00020',
'dG_mr0_ml1_ev0p00005_sO_c50_mr0_ml1_ev0p00005_hinsen'  ,
'dG_mr1_ml0_ev0p00020_sO_c50_mr1_ml0_ev0p00020_hinsen',
'dG_mr0_ml1_ev0p00008_sO_c50_mr0_ml1_ev0p00008'  ,
'dG_mr1_ml0_ev0p00030_sO_c50_mr1_ml0_ev0p00030',
'dG_mr0_ml1_ev0p00008_sO_c50_mr0_ml1_ev0p00008_hinsen'  ,
'dG_mr1_ml0_ev0p00030_sO_c50_mr1_ml0_ev0p00030_hinsen',
'dG_mr0_ml1_ev0p00010_sO_c50_mr0_ml1_ev0p00010'  ,
'dG_mr1_ml0_ev0p00035_sO_c50_mr1_ml0_ev0p00035',
'dG_mr0_ml1_ev0p00010_sO_c50_mr0_ml1_ev0p00010_hinsen'  ,
'dG_mr1_ml0_ev0p00035_sO_c50_mr1_ml0_ev0p00035_hinsen',
'dG_mr0_ml1_ev0p00012_sO_c50_mr0_ml1_ev0p00012'  ,
'dG_mr1_ml0_ev0p00050_sO_c50_mr1_ml0_ev0p00050',
'dG_mr0_ml1_ev0p00012_sO_c50_mr0_ml1_ev0p00012_hinsen'  ,
'dG_mr1_ml0_ev0p00050_sO_c50_mr1_ml0_ev0p00050_hinsen',
'dG_mr0_ml1_ev0p00014_sO_c50_mr0_ml1_ev0p00014'  ,
'dG_mr1_ml0_ev0p00080_sO_c50_mr1_ml0_ev0p00080',
'dG_mr0_ml1_ev0p00014_sO_c50_mr0_ml1_ev0p00014_hinsen'  ,
'dG_mr1_ml0_ev0p00080_sO_c50_mr1_ml0_ev0p00080_hinsen',
'dG_mr0_ml1_ev0p00015_sO_c50_mr0_ml1_ev0p00015'  ,
'dG_mr1_ml0_ev0p70000_sO_c50_mr1_ml0_ev0p70000',
'dG_mr0_ml1_ev0p00015_sO_c50_mr0_ml1_ev0p00015_hinsen'  ,
'dG_mr1_ml0_ev1p00000_sO_c50_mr1_ml0_ev1p00000',
'dG_mr0_ml1_ev0p00016_sO_c50_mr0_ml1_ev0p00016'  ,
'dG_mr1_ml0_ev1p30000_sO_c50_mr1_ml0_ev1p30000',
'dG_mr0_ml1_ev0p00016_sO_c50_mr0_ml1_ev0p00016_hinsen'  ,
'dG_mr1_ml0_ev1p70000_sO_c50_mr1_ml0_ev1p70000',
'dG_mr0_ml1_ev0p00017_sO_c50_mr0_ml1_ev0p00017'  ,
'dG_mr1_ml0_ev2p00000_sO_c50_mr1_ml0_ev2p00000',
'dG_mr0_ml1_ev0p00017_sO_c50_mr0_ml1_ev0p00017_hinsen'  ,
'dG_mr1_ml0_ev2p50000_sO_c50_mr1_ml0_ev2p50000',
'dG_mr0_ml1_ev0p00018_sO_c50_mr0_ml1_ev0p00018'  ,
'dG_mr1_ml1_ev0p00005_sO_c50_mr1_ml1_ev0p00005',
'dG_mr0_ml1_ev0p00018_sO_c50_mr0_ml1_ev0p00018_hinsen'  ,
'dG_mr1_ml1_ev0p00005_sO_c50_mr1_ml1_ev0p00005_hinsen',
'dG_mr0_ml1_ev0p00020_sO_c50_mr0_ml1_ev0p00020'  ,
'dG_mr1_ml1_ev0p00008_sO_c50_mr1_ml1_ev0p00008',
'dG_mr0_ml1_ev0p00020_sO_c50_mr0_ml1_ev0p00020_hinsen'  ,
'dG_mr1_ml1_ev0p00008_sO_c50_mr1_ml1_ev0p00008_hinsen',
'dG_mr0_ml1_ev0p00030_sO_c50_mr0_ml1_ev0p00030'  ,
'dG_mr1_ml1_ev0p00010_sO_c50_mr1_ml1_ev0p00010',
'dG_mr0_ml1_ev0p00030_sO_c50_mr0_ml1_ev0p00030_hinsen'  ,
'dG_mr1_ml1_ev0p00010_sO_c50_mr1_ml1_ev0p00010_hinsen',
'dG_mr0_ml1_ev0p00035_sO_c50_mr0_ml1_ev0p00035'  ,
'dG_mr1_ml1_ev0p00012_sO_c50_mr1_ml1_ev0p00012',
'dG_mr0_ml1_ev0p00035_sO_c50_mr0_ml1_ev0p00035_hinsen'  ,
'dG_mr1_ml1_ev0p00012_sO_c50_mr1_ml1_ev0p00012_hinsen',
'dG_mr0_ml1_ev0p00050_sO_c50_mr0_ml1_ev0p00050'  ,
'dG_mr1_ml1_ev0p00014_sO_c50_mr1_ml1_ev0p00014',
'dG_mr0_ml1_ev0p00050_sO_c50_mr0_ml1_ev0p00050_hinsen'  ,
'dG_mr1_ml1_ev0p00014_sO_c50_mr1_ml1_ev0p00014_hinsen',
'dG_mr0_ml1_ev0p00080_sO_c50_mr0_ml1_ev0p00080'  ,
'dG_mr1_ml1_ev0p00015_sO_c50_mr1_ml1_ev0p00015',
'dG_mr0_ml1_ev0p00080_sO_c50_mr0_ml1_ev0p00080_hinsen'  ,
'dG_mr1_ml1_ev0p00015_sO_c50_mr1_ml1_ev0p00015_hinsen',
'dG_mr0_ml1_ev0p70000_sO_c50_mr0_ml1_ev0p70000'  ,
'dG_mr1_ml1_ev0p00016_sO_c50_mr1_ml1_ev0p00016',
'dG_mr0_ml1_ev1p00000_sO_c50_mr0_ml1_ev1p00000'  ,
'dG_mr1_ml1_ev0p00016_sO_c50_mr1_ml1_ev0p00016_hinsen',
'dG_mr0_ml1_ev1p30000_sO_c50_mr0_ml1_ev1p30000'  ,
'dG_mr1_ml1_ev0p00017_sO_c50_mr1_ml1_ev0p00017',
'dG_mr0_ml1_ev1p70000_sO_c50_mr0_ml1_ev1p70000'  ,
'dG_mr1_ml1_ev0p00017_sO_c50_mr1_ml1_ev0p00017_hinsen',
'dG_mr0_ml1_ev2p00000_sO_c50_mr0_ml1_ev2p00000'  ,
'dG_mr1_ml1_ev0p00018_sO_c50_mr1_ml1_ev0p00018',
'dG_mr0_ml1_ev2p50000_sO_c50_mr0_ml1_ev2p50000'  ,
'dG_mr1_ml1_ev0p00018_sO_c50_mr1_ml1_ev0p00018_hinsen',
'dG_mr0_ml3_ev0p00015_sO_c50_mr0_ml3_ev0p00015hinsen'  ,
'dG_mr1_ml1_ev0p00020_sO_c50_mr1_ml1_ev0p00020',
'dG_mr0_ml3_ev0p00016_sO_c50_mr0_ml3_ev0p00016hinsen'  ,
'dG_mr1_ml1_ev0p00020_sO_c50_mr1_ml1_ev0p00020_hinsen',
'dG_mr0_ml3_ev0p00017_sO_c50_mr0_ml3_ev0p00017hinsen'  ,
'dG_mr1_ml1_ev0p00030_sO_c50_mr1_ml1_ev0p00030',
'dG_mr0_ml3_ev0p00018_sO_c50_mr0_ml3_ev0p00018hinsen'  ,
'dG_mr1_ml1_ev0p00030_sO_c50_mr1_ml1_ev0p00030_hinsen',
'dG_mr0_ml3_ev0p00020_sO_c50_mr0_ml3_ev0p00020hinsen'  ,
'dG_mr1_ml1_ev0p00035_sO_c50_mr1_ml1_ev0p00035',
'dG_mr10_ml0_ev0p00005_sO_c50_mr10_ml0_ev0p00005hinsen'  ,
'dG_mr1_ml1_ev0p00035_sO_c50_mr1_ml1_ev0p00035_hinsen',
'dG_mr10_ml0_ev0p00008_sO_c50_mr10_ml0_ev0p00008hinsen'  ,
'dG_mr1_ml1_ev0p00050_sO_c50_mr1_ml1_ev0p00050',
'dG_mr10_ml0_ev0p00010_sO_c50_mr10_ml0_ev0p00010hinsen'  ,
'dG_mr1_ml1_ev0p00050_sO_c50_mr1_ml1_ev0p00050_hinsen',
'dG_mr10_ml0_ev0p00012_sO_c50_mr10_ml0_ev0p00012hinsen'  ,
'dG_mr1_ml1_ev0p00080_sO_c50_mr1_ml1_ev0p00080',
'dG_mr10_ml0_ev0p00014_sO_c50_mr10_ml0_ev0p00014hinsen'  ,
'dG_mr1_ml1_ev0p00080_sO_c50_mr1_ml1_ev0p00080_hinsen',
'dG_mr10_ml0_ev0p00015_sO_c50_mr10_ml0_ev0p00015hinsen'  ,
'dG_mr1_ml1_ev0p70000_sO_c50_mr1_ml1_ev0p70000',
'dG_mr10_ml0_ev0p00016_sO_c50_mr10_ml0_ev0p00016hinsen'  ,
'dG_mr1_ml1_ev1p00000_sO_c50_mr1_ml1_ev1p00000',
'dG_mr10_ml0_ev0p00017_sO_c50_mr10_ml0_ev0p00017hinsen'  ,
'dG_mr1_ml1_ev1p30000_sO_c50_mr1_ml1_ev1p30000',
'dG_mr10_ml0_ev0p00018_sO_c50_mr10_ml0_ev0p00018hinsen'  ,
'dG_mr1_ml1_ev1p70000_sO_c50_mr1_ml1_ev1p70000',
'dG_mr10_ml0_ev0p00020_sO_c50_mr10_ml0_ev0p00020hinsen'  ,
'dG_mr1_ml1_ev2p00000_sO_c50_mr1_ml1_ev2p00000',
'dG_mr10_ml0_ev0p00030_sO_c50_mr10_ml0_ev0p00030hinsen'  ,
'dG_mr1_ml1_ev2p50000_sO_c50_mr1_ml1_ev2p50000',
'dG_mr10_ml0_ev0p00035_sO_c50_mr10_ml0_ev0p00035hinsen'  ,
'dG_mr3_ml0_ev0p00015_sO_c50_mr3_ml0_ev0p00015hinsen',
'dG_mr10_ml0_ev0p00050_sO_c50_mr10_ml0_ev0p00050hinsen'  ,
'dG_mr3_ml0_ev0p00016_sO_c50_mr3_ml0_ev0p00016hinsen',
'dG_mr10_ml0_ev0p00080_sO_c50_mr10_ml0_ev0p00080hinsen'  ,
'dG_mr3_ml0_ev0p00017_sO_c50_mr3_ml0_ev0p00017hinsen',
'dG_mr10_ml0_ev0p70000_sO_c50_mr10_ml0_ev0p70000'  ,
'dG_mr3_ml0_ev0p00018_sO_c50_mr3_ml0_ev0p00018hinsen',
'dG_mr10_ml0_ev1p00000_sO_c50_mr10_ml0_ev1p00000'  ,
'dG_mr3_ml0_ev0p00020_sO_c50_mr3_ml0_ev0p00020hinsen',
'dG_mr10_ml0_ev1p30000_sO_c50_mr10_ml0_ev1p30000'  ,
'dG_mr3_ml3_ev0p00015_sO_c50_mr3_ml3_ev0p00015hinsen',
'dG_mr10_ml0_ev1p70000_sO_c50_mr10_ml0_ev1p70000'  ,
'dG_mr3_ml3_ev0p00016_sO_c50_mr3_ml3_ev0p00016hinsen',
'dG_mr10_ml0_ev2p00000_sO_c50_mr10_ml0_ev2p00000'  ,
'dG_mr3_ml3_ev0p00017_sO_c50_mr3_ml3_ev0p00017hinsen',
'dG_mr10_ml0_ev2p50000_sO_c50_mr10_ml0_ev2p50000'  ,
'dG_mr3_ml3_ev0p00018_sO_c50_mr3_ml3_ev0p00018hinsen',
'dG_mr10_ml10_ev0p00005_sO_c50_mr10_ml10_ev0p00005hinsen'  ,
'dG_mr3_ml3_ev0p00020_sO_c50_mr3_ml3_ev0p00020hinsen',
'dG_mr10_ml10_ev0p00008_sO_c50_mr10_ml10_ev0p00008hinsen'  ,
'dG_mr10_ml10_ev0p00010_sO_c50_mr10_ml10_ev0p00010hinsen'  ,
'dG_mr10_ml10_ev0p00012_sO_c50_mr10_ml10_ev0p00012hinsen'  ,
'dG_mr10_ml10_ev0p00014_sO_c50_mr10_ml10_ev0p00014hinsen'  ,
'dG_mr10_ml10_ev0p00015_sO_c50_mr10_ml10_ev0p00015hinsen'  ,
'dG_mr10_ml10_ev0p00016_sO_c50_mr10_ml10_ev0p00016hinsen'  ,
'dG_mr10_ml10_ev0p00017_sO_c50_mr10_ml10_ev0p00017hinsen'  ,
'dG_mr10_ml10_ev0p00018_sO_c50_mr10_ml10_ev0p00018hinsen'  ,
'dG_mr10_ml10_ev0p00020_sO_c50_mr10_ml10_ev0p00020hinsen',
'dG_mr3_ml0_ev0p00021_sO_c50_mr3_ml0_ev0p00021hinsen',
'dG_mr3_ml0_ev0p00022_sO_c50_mr3_ml0_ev0p00022hinsen',
'dG_mr3_ml0_ev0p00023_sO_c50_mr3_ml0_ev0p00023hinsen',
'dG_mr3_ml0_ev0p00024_sO_c50_mr3_ml0_ev0p00024hinsen',
'dG_mr3_ml0_ev0p00025_sO_c50_mr3_ml0_ev0p00025hinsen',
'dG_mr5_ml0_ev0p00025_sO_c50_mr5_ml0_ev0p00025hinsen',
'dG_mr5_ml0_ev0p00027_sO_c50_mr5_ml0_ev0p00027hinsen',
'dG_mr5_ml0_ev0p00029_sO_c50_mr5_ml0_ev0p00029hinsen',
'dG_mr5_ml0_ev0p00030_sO_c50_mr5_ml0_ev0p00030hinsen',
'dG_mr5_ml0_ev0p00035_sO_c50_mr5_ml0_ev0p00035hinsen'
'dG_mr5_ml0_ev0p00038_sO_c50_mr5_ml0_ev0p00038hinsen'
]


evalData(protList, benchmarkList,path, filename)

def getProtein(dataframe, protein):
    return dataframe.loc[dataframe['protein']== protein]

def getHinsen(dataframe, brownian):
    return dataframe.loc[dataframe['brownian']== brownian]

def getEVScale(dataframe, evScale):
    return dataframe.loc[dataframe['evScale']== evScale]

def getModesRec(dataframe, numModes):
    return dataframe.loc[dataframe['numModesRec']== numModes]

def getModesLig(dataframe, numModes):
    return dataframe.loc[dataframe['numModesLig']== numModes]






def readDataframe(filename):
    return pd.read_csv(filename)


zeroframe = pf.loc[(pf['numModesRec'] == 0) &( pf['numModesLig'] == 0) & (pf['evScale'] ==1) & (pf['brownian'] == 0)]

def createDiffFrame(frame, refFrame):
    pf = frame
    #pf['sum'] = pf['one_star_50']+pf['two_star_50']+pf['three_star_50']
    pf['sum'] = pf['one_star_50']+2*pf['two_star_50']+3*pf['three_star_50']
    zeroframe = refFrame
    pf['diff'] = 0
    pf['diff_one'] = 0
    pf['diff_two']  = 0
    pf['diff_three'] = 0

    for idx in pf.index:
        protein = pf.at[idx, 'protein']
        totdiff = pf.at[idx, 'sum'] - zeroframe.loc[zeroframe['protein'] == protein]['sum']
        pf._set_value(idx, 'diff', totdiff)

        onediff = pf.at[idx, 'one_star_50'] - zeroframe.loc[zeroframe['protein'] == protein]['one_star_50']
        pf._set_value(idx, 'diff_one', onediff)
        twodiff = pf.at[idx, 'two_star_50'] - zeroframe.loc[zeroframe['protein'] == protein]['two_star_50']
        pf._set_value(idx, 'diff_two', twodiff)
        threediff = pf.at[idx, 'three_star_50'] - zeroframe.loc[zeroframe['protein'] == protein]['three_star_50']
        pf._set_value(idx, 'diff_three', threediff)
    return pf


def binProteins(dataframe):
    list = []
    pf = dataframe
    brownian = pf['brownian'].unique()
    for i in brownian:
        scales = pf.loc[pf['brownian'] == i]['evScale'].unique()
        modes = pf.loc[(pf['brownian'] == i) ]['numModesRec'].unique() #& (pf['evScale'] == scale)
        for scale in scales:
            for mode in modes:
                df = pf.loc[(pf['numModesRec'] == mode) & (pf['numModesLig'] == mode) & (pf['evScale'] == scale) & (
                        pf['brownian'] == i)]
                one =df['one_star_50'].sum()
                two = df['two_star_50'].sum()
                three = df['three_star_50'].sum()

                row = {'numModes': mode, 'scale': scale, 'brownian': i, 'one': one, 'two': two, 'three': three,
                       'sum': one + two + three}
                list.append(row)
    return    pd.DataFrame(list)



def plotProteinBins(dataframe):
    plt.clf()
    benchmarkComp = dataframe
    fig, axes = plt.subplots(4, 2)
    brownian = benchmarkComp['brownian'].unique()

    for i in brownian:
        scales = benchmarkComp.loc[benchmarkComp['brownian'] == i]['scale'].unique()
        for scale in scales:
            sub = benchmarkComp.loc[(benchmarkComp['scale'] == scale) & (benchmarkComp['brownian'] == i)].sort_values('numModes')
            if i == 1:
                add = 0.5
            else:
                add = 0
            axes[0,i].plot(sub['numModes'], sub['one'], label = "s {} b {} ".format(scale, i))
            axes[1,i].plot(sub['numModes'], sub['two'], label="s {} b {} ".format(scale, i))
            axes[2,i].plot(sub['numModes'], sub['three'], label="s {} b {} ".format(scale, i))
            axes[3, i].plot(sub['numModes'], sub['sum'], label="s {} b {} ".format(scale, i))
        axes[0,i].legend()
        axes[1,i].legend()
        axes[2,i].legend()
        axes[3,i].legend()


def createSum(df):
    df['sum'] = df['one_star_50']+ 2*df['two_star_50']+ 3*df['three_star_50']
    return df


# pf.loc[pf['protein'] != '7CEI'].nsmallest(200,'diff')[['protein', 'numModesRec','numModesLig','evScale','brownian','diff']]
# pf.nsmallest(1000,'diff')['protein'].value_counts()
# pf.nlargest(1000,'diff')['protein'].value_counts()
# pf.loc[df['brownian'] == 1].nsmallest(100,'diff')['numModesRec'].plot.hist(bins = 100)
#
# sub = pf.sort_values(['three_star_50','two_star_50','one_star_50'], ascending = False)[:10][
#     ['one_star_50','two_star_50','three_star_50']]


# fileList = ["/home/glenn/Documents/Masterarbeit/analysis/largerun/benchmark5_best1","/home/glenn/Documents/Masterarbeit/analysis/largerun/benchmark5_worst1"]
# for filename in fileList:
#     for brownian in [0,1]:
#         #filename = "/home/glenn/Documents/Masterarbeit/analysis/largerun/benchmark5_best"
#         df = readDataframe(filename)
#         df = createDiffFrame(df)
#         binned = binProteins(df)
#         plotProteinBins(binned)
#         brownian = 0
#         n = 50
#          print( "\n" ,df.nlargest(100,'diff')['brownian'].value_counts())
#          print( "\n" , n, "largest ", " protein  \n" ,df.loc[df['brownian'] == brownian].nlargest(n,'diff')['protein'].value_counts())
#          print( "\n" , n, "largest ", " modesrec \n" ,df.loc[df['brownian'] == brownian].nlargest(n,'diff')['numModesRec'].value_counts())
#          print( "\n" , n, "largest ", " modeslig \n" ,df.loc[df['brownian'] == brownian].nlargest(n,'diff')['numModesLig'].value_counts())
#          print( "\n" , n, "largest ", " scale    \n" ,df.loc[df['brownian'] == brownian].nlargest(n,'diff')['evScale'].value_counts())
#          print( "\n" , n, "smallest ", " protein  \n" , df.loc[df['brownian'] == brownian].nsmallest(n,'diff')['protein'].value_counts())
#          print( "\n" , n, "smallest ", " modesrec \n" , df.loc[df['brownian'] == brownian].nsmallest(n,'diff')['numModesRec'].value_counts())
#          print( "\n" , n, "smallest ", " modeslig \n" , df.loc[df['brownian'] == brownian].nsmallest(n,'diff')['numModesLig'].value_counts())
#          print( "\n" , n, "smallest ", " scale    \n" , df.loc[df['brownian'] == brownian].nsmallest(n,'diff')['evScale'].value_counts())
#
# plt.clf()
# brownian =1
# col = "evScale"
# subset = df.loc[df['brownian'] == brownian].nlargest(n,'diff'
# )[col].value_counts()
# plt.scatter(subset.index, subset.values)
# plt.xlabel(col)
# plt.ylabel("Counts in top 50 of sum")
# plt.savefig("/home/glenn/Documents/Masterarbeit/analysis/benchhmark5_best/largest_value_counts_brownian_{}_col_{}.png".format(brownian,col))

# subset = df.loc[df['brownian'] == brownian].nlargest(n,'diff')['protein'].value_counts())
# plt.scatter(subset.index, subset.values)
# plt.show()

fext['pdb'] = '-receptor-for-docking'
protList = readStringList("/home/glenn/work/benchmark5_attract/protlist",0)
dele = ['1QFW','1BJ1', '1FSK', '1IQD', '1K4C', '1KXQ', '1NCA', '1NSN', '1QFW', '2HMI', '2JEL', '9QFW','1DE4','4FQI','4GAM','4GXU','1EXB','4FQI' ,'1EER', '4H03']
for d in dele:
    protlist.remove(d)

# dm = loadRmsd(protList,['dG_mr5_ml5_ev1p0_sO_c50_mr5_ml5_ev1'],'/home/glenn/work/benchmark5_attract' )
#
# dm = loadRmsd(protList,['dG_mr10_ml10_ev1p0_sO_c50_mr10_ml10_ev1p0'],'/home/glenn/work/benchmark5_attract' )
# d = loadRmsd(protList,['dG_mr0_ml0_ev1p0_sO_c50_mr0_ml0_ev1'],'/home/glenn/work/benchmark5_attract' )
# proteins= list(set(list(d)) & set(list(dm)))
#
# dtest= dm
# rmsdref = d
#
# def listToDict(list):
#     return {i: x for i, x in enumerate(list)}
#
# def getLargestRMSDIdx(searchlist, topnumber):
#     return np.argsort(-searchlist)[:topnumber]
#
# def getIndicesOfSimilarEntries(list, searchobj):
#     return np.where(np.isclose(list, searchobj, 1e-03, 1e-03))[0]
#
# def getBestMatch(list, indices, referenceobj):
#     bestmatch = 100
#     bestidx = 0
#     for idx in indices:
#         if abs(list[idx] - referenceobj) < bestmatch:
#             bestmatch = ref[idx]
#             bestidx = idx
#     return bestmatch, bestidx
#
#
# count = 0
# num = 20
# minnum = 10
# resList =[]
# for protein in proteins:
#     tlist = dtest[protein][:num]
#     lagestRmsdIdx = getLargestRMSDIdx(tlist, minnum)
#     ref = rmsdref[protein][:num]
#     for index in lagestRmsdIdx:
#         tobj = tlist[index]
#         indices = getIndicesOfSimilarEntries(ref, tobj)
#         if len(indices) > 0:
#             count += 1
#             match, matchIndex = getBestMatch(ref, indices, tobj)
#             #print match , tobj
#             #if ( index-matchIndex) < 0:
#             print protein,"best match" ,match ,"val",tobj, "bidx",matchIndex ,"idx ",index, index-matchIndex
#             result = {"protein":protein,"value": tobj, "valIdx": index, "match": match, "mIdx":matchIndex, "idff":index-matchIndex}
#             resList.append(result)
# resF = pd.DataFrame(resList)
#
# print count
# print resF.nsmallest(50,'idff')
# listLargestDiff = resF.nsmallest(50,'idff')['protein'].values
# print resF.nsmallest(50,'idff')['protein'].value_counts().nlargest(10)
# listCounts = resF.nsmallest(50,'idff')['protein'].value_counts().nlargest(10).index.tolist()




#saveRmsdEnergy(listCounts, ['dG_mr0_ml0_ev1p0_sO_c50_mr0_ml0_ev1','dG_mr5_ml5_ev1p0_sO_c50_mr5_ml5_ev1'],"/home/glenn/work/benchmark5_attract",'/home/glenn/Documents/Masterarbeit/analysis/largerun/rmsdplot/dG_mr0_ml0_ev1p0_sO_c50_mr0_ml0_ev1_VS_dG_mr5_ml5_ev1p0_sO_c50_mr5_ml5_ev___MOSTCOUNTS')
# for prot in listCounts:
#     plt.waitforbuttonpress()
#     plt.clf()
#     plotRmsdEnergy([prot], ['dG_mr10_ml10_ev1p0_sO_c50_mr10_ml10_ev1p0','dG_mr5_ml5_ev1p0_sO_c50_mr5_ml5_ev1','dG_mr0_ml0_ev1p0_sO_c50_mr0_ml0_ev1'], "/home/glenn/work/benchmark5_attract")
#     plt.xlim([1,80])
#     plt.ylim([-30, 5])
#     plt.show()
#
# for prot in listLargestDiff:
#     plt.waitforbuttonpress()
#     plt.clf()
#     plotRmsdEnergy([prot], ['dG_mr10_ml10_ev1p0_sO_c50_mr10_ml10_ev1p0','dG_mr5_ml5_ev1p0_sO_c50_mr5_ml5_ev1','dG_mr0_ml0_ev1p0_sO_c50_mr0_ml0_ev1'], "/home/glenn/work/benchmark5_attract")
#     #plt.xlim([1,80])
#     plt.ylim([-30, 5])
#     #plt.savefig('/home/glenn/Documents/Masterarbeit/analysis/rmsdworsening/{}.png'.format(prot))
#     plt.show()

#diff bound and omodes
# df_omodes['d_bound'] = 0
# sum1 = 0; sum2 = 0; sum3 = 0
# for i in df_omodes.loc[df_omodes['evScale'] == 1.0].index:
#     prot = df_omodes.at[i,'protein']
#     rowb = df_bound.loc[df_bound['protein']== prot]
#     b1 = rowb['one_star_50']
#     b2 = rowb['two_star_50']
#     b3 = rowb['three_star_50']
#     u1 = df_omodes.at[i,'one_star_50']
#     u2 = df_omodes.at[i,'two_star_50']
#     u3 = df_omodes.at[i,'three_star_50']
#     d1 = u1 - b1
#     d2 = u2 - b2
#     d3 = u3 - b3
#     df_omodes.at[i,'d_bound'] = u1-b1+2*(u2-b2)+3*(u3-b3)
#     print "bound " ,b1.values[0], b2.values[0], b3.values[0], "unbound" , u1,u2,u3
#     if b1.values[0] < u1 and b2.values[0] < u2 and b3.values[0] < u3:
#         print prot , "unbound besser"
#     sum1 += d1.values[0]
#     sum2 += d2.values[0]
#     sum3 += d3.values[0]
# print sum1, sum2, sum3
# df_omodes.nlargest(10,'d_bound')['protein']


#good:
#1AK4
#1AZS
#bad:
#1AY7
#1BKD
# rmsdListHigh=['1HCF','1QFW','1WDW','1PXV','1DFJ','1XD3','4CPA','1EWY','1MAH','2OT3','1KXP','4DN4','1HCF','3H2V']
# goodProt=['1AK4', "3PC8"]
# listtest = ['3BX7']
# for prot in listtest:
#     #plt.waitforbuttonpress()
#     plt.clf()
#
#     plotRmsdEnergy([prot], ['dG_mr10_ml1_ev1p0_sO_c50_mr10_ml1_ev1p0','dG_mr1_ml0_ev1p0_sO_c50_mr1_ml0_ev1p0',
#                             'dG_mr0_ml0_ev1p0_sO_c50_mr0_ml0_ev1p0'], "/home/glenn/work/benchmark5_attract")
#     plt.xlim([1,100])
#     plt.ylim([-30, 5])
#     plt.ylabel("Energy")
#     plt.title("Protein {}".format(prot))
#     plt.xlabel("RMSD")
#     #plt.savefig('/home/glenn/Documents/Masterarbeit/analysis/rmsdworsening/best_{}.png'.format(prot))
#     #plt.show()
#
#
# di = getEnergyDelta(protList, ['dG_mr10_ml1_ev1p0_sO_c50_mr10_ml1_ev1p0',
#                             'dG_mr0_ml0_ev1p0_sO_c50_mr0_ml0_ev1p0'], "/home/glenn/work/benchmark5_attract", 100)