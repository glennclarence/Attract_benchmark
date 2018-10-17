import os
#from  load_pdbs import ProteinEsemble
import numpy as np
import pandas as pd

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
                data[i] = np.fromstring(line, dtype=type, sep=sep)[colIdx]
                #data[i] = type(line.split()[colIdx])
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
            energies = 0
        return np.asarray( energies )
    except:
        return -1

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


bm="dG_mr5_ml1_ev0p1_sO_c50_mr5_ml1_ev0p1_hinsen"
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

benchmarkList = ['dG_mr5_ml0_ev1p0_sO_c50_mr5_ml0_ev1p0']
protList_best = [
'1AK4',
'1IJK',
'1K4C',
'1K74',
'1KXQ',
'1RV6',
'1SYX',
'1XQS',
'3MXW',
'3PC8',
'4G6M'
]

protList_worst =  [
'1AVX',
'1B6C',
'1BUH',
'1GCQ',
'1I9R',
'1PPE',
'1ZHI',
'2FD6',
'2OUL',
'2SNI',
'7CEI'
]
basePath ='/home/glenn/work/benchmark5_worst'

def readStringList(filename, colIdx):
    with open(filename) as f:
        lines = f.readlines()
        list = []
        for line in lines:
            list.append(line.split()[colIdx])
        return  list

benchmarkList = readStringList('/home/glenn/work/benchmarks_worst',0)


def isBrownian(bmname):
    return
#protList = ['3MXW']
#benchmarkList = ['dG_mr5_ml0_ev0p05_sO_c50_mr5_ml0_ev0p05_hinsen']
# dataList=[]
# for protein in protList:
#     for bm in benchmarkList:
#         bmPath = Path(getPath(basePath=basePath, benchmarkName=bm, proteinName= protein))
#         evScale = getModeEVScale(bm)
#         numModesRec = getNumModesFromString(bm, 'r')
#         numModesLig = getNumModesFromString(bm, 'l')
#         baseFile = protein+fext['pdb']
#
#         energyFile = baseFile + fext['energy']
#         if(bmPath.isfile(energyFile)):
#             energies = getEnergyfromFile(bmPath.filename(energyFile))
#         rmsdFile = baseFile + fext['rmsd']
#         if (bmPath.isfile(rmsdFile)):
#             rmsd = np.asarray(readRMSD(bmPath.filename(rmsdFile)))
#         irmsdFile = baseFile + fext['irmsd']
#         if (bmPath.isfile(irmsdFile)):
#             irmsd = np.asarray(readRMSD(bmPath.filename(irmsdFile)))
#         fnatFile = baseFile + fext['fnat']
#         if (bmPath.isfile(fnatFile)):
#             fnat = np.asarray(readFnat(bmPath.filename(fnatFile)))
#         if type(fnat) == int or type(irmsd) == int or type(rmsd) == int :
#             print "file not available for ", protein
#             continue
#         capriList = getCapriList(irmsd, fnat)
#         one_star, two_star, three_star = capriCount(capriList)
#         print bm.find('hinsen')
#         if bm.find('hinsen') != -1:
#             brownian = 1
#         else:
#             brownian = 0
#         row = {'protein':protein,'numModesRec':numModesRec,'numModesLig':numModesLig,'evScale':evScale, 'one_star_50':one_star,'two_star_50':two_star,'three_star_50':three_star, 'brownian': brownian}
#         dataList.append(row)
#         print protein,numModesLig, numModesRec, evScale,one_star, two_star, three_star, brownian
# frame = pd.DataFrame(dataList)
# frame.to_csv("worst_csv")
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



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

brownian = 0
scale = 1.0
frame = pd.read_csv('/home/glenn/Documents/Masterarbeit/analysis/bestProteins/evaluation/data_benchmark_best_hinsen_noHinsen')
protein = '1AK4'
plt.clf()
#plt.gca().set_prop_cycle(None)
pf = getProtein(frame, protein)
scales = pf['evScale'].loc[pf['brownian']==brownian].unique()
print "scales", scales
xticks= []

for i,scale in enumerate(scales):

    one = pf.loc[(pf['numModesRec'] == pf['numModesLig']) & ( pf['evScale'] == scale) & (pf['brownian']== brownian)]['one_star_50']
    two = pf.loc[(pf['numModesRec'] == pf['numModesLig']) & ( pf['evScale'] == scale) & (pf['brownian']== brownian)]['two_star_50']
    three = pf.loc[(pf['numModesRec'] == pf['numModesLig']) & ( pf['evScale'] == scale) & (pf['brownian']== brownian)]['three_star_50']

    segSize = len(one)

    if segSize == 0:
        continue
    seg = np.arange(segSize)
    print scale, i, segSize
    numModes = pf.loc[(pf['numModesRec'] == pf['numModesLig']) & ( pf['evScale'] == scale) & (pf['brownian']== brownian)]['numModesRec']
    ticks = ['%s\n%s' % (k,scale) for k in numModes.values.tolist()]
    xticks = xticks +ticks
    seg += i*segSize

    plt.bar(seg,three, width = 1,alpha = 0.5,color='red')
    plt.bar(seg,two, width = 1, bottom = three ,alpha = 0.5,color='blue')
    plt.bar(seg,one, width = 1, bottom = two,alpha = 0.5,color='magenta')
plt.xticks(np.arange(len(xticks)),xticks)
    # color=['red', 'green', 'blue', 'cyan', 'magenta'],


protList = protList_best
filename= '/home/glenn/Documents/Masterarbeit/analysis/bestProteins/evaluation/data_benchmark_best_hinsen_noHinsen'
frame = pd.read_csv(filename)

list = []
pf = frame

brownian =[0,1]
for i in brownian:
    scales = pf.loc[pf['brownian'] == i]['evScale'].unique()

    modes = pf.loc[(pf['brownian'] == i) ]['numModesRec'].unique() #& (pf['evScale'] == scale)
    print scales, modes
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
benchmarkComp = pd.DataFrame(list)
benchmarkComp.loc[benchmarkComp['brownian'] == 0]['scale'].unique()

plt.clf()
fig, axes = plt.subplots(4,2)
for i in brownian:
    scales = pf.loc[pf['brownian'] == i]['evScale'].unique()

    modes = pf.loc[(pf['brownian'] == i)]['numModesRec'].unique()  # & (pf['evScale'] == scale)
    for scale in scales:
        sub = benchmarkComp.loc[(benchmarkComp['scale'] == scale) & (benchmarkComp['brownian'] == i)].sort_values('numModes')
        if brownian == 1:
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



protList = protList_worst
filename= '/home/glenn/Documents/Masterarbeit/analysis/worst_csv'
frame = pd.read_csv(filename)

list = []
pf = frame
pf['sum'] = pf['one_star_50']+pf['two_star_50']+pf['three_star_50']

zeroframe = pf.loc[(pf['numModesRec'] == 0) &( pf['numModesLig'] == 0) & (pf['evScale'] ==1) & (pf['brownian'] == 0)]

brownian =[0,1]
for i in brownian:
    scales = pf.loc[pf['brownian'] == i]['evScale'].unique()

    modes = pf.loc[(pf['brownian'] == i) ]['numModesRec'].unique() #& (pf['evScale'] == scale)
    print scales, modes
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
benchmarkComp = pd.DataFrame(list)
benchmarkComp.loc[benchmarkComp['brownian'] == 0]['scale'].unique()

plt.clf()
fig, axes = plt.subplots(4,2)
for i in brownian:
    scales = pf.loc[pf['brownian'] == i]['evScale'].unique()

    modes = pf.loc[(pf['brownian'] == i)]['numModesRec'].unique()  # & (pf['evScale'] == scale)
    for scale in scales:
        sub = benchmarkComp.loc[(benchmarkComp['scale'] == scale) & (benchmarkComp['brownian'] == i)].sort_values('numModes')
        if brownian == 1:
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





sub = pf.sort_values(['three_star_50','two_star_50','one_star_50'], ascending = False)[:10][
    ['one_star_50','two_star_50','three_star_50']]
