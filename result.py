import os
#from  load_pdbs import ProteinEsemble
import numpy as np
import matplotlib.pyplot as plt
from eval_class import ResultClass
from collections import OrderedDict

# idx_mode = 0
# idx_rmsd = 1
# idx_energy = 2
# idx_pos = 3
# idx_min = 4
# idx_mean = 5
# idx_mean_10 = 6
# idx_mean_50 = 7
# idx_mean_sort_10 = 8
# idx_mean_sort_50 = 9
# idx_amp = 10
# idx_ev = 11
# tot_num = 12

dict_indices = {}
dict_indices['mode'] = 0
dict_indices['rmsd'] = 1
dict_indices['energy'] = 2
dict_indices['pos'] = 3
dict_indices['min'] = 4
dict_indices['mean'] = 5
dict_indices['mean_10'] = 6
dict_indices['mean_50'] = 7
dict_indices['mean_sort_10'] = 8
dict_indices['mean_sort_50'] = 9
dict_indices['amp'] = 10
dict_indices['ev'] = 11


#filename_scoring = "/home/glenn/cluster/benchmark_attract_test/1AVX/benchmark_ORI_scorig_50cut_5modes_2/1AVX-receptor-for-docking-sorted-dr.dat"
#energies =  getEnergyfromFile( filename_scoring)
#print energies

a5BM = ResultClass( dict_indices )

path_folder = "/home/glenn/Documents/Masterarbeit/test_bencheval/"
path_folder= "/home/glenn/cluster/benchmark5_attract"
BM5m = "benchmark_GPU_scorig_50cut_5modes"
BM0m = "benchmark_GPU_scorig_50cut_0modes"
BM3m = "benchmark_GPU_scorig_50cut_3modes"


BMLoad = [
'benchmark_GPU_scGPU_50cut_0modes',
'benchmark_GPU_scGPU_50cut_3modes',
'benchmark_GPU_scGPU_50cut_3modes_sc2EV',
'benchmark_GPU_scGPU_50cut_3modes_scnoEV',
'benchmark_GPU_scGPU_50cut_5modes',
'benchmark_GPU_scGPU_50cut_5modes_sc2EV',
'benchmark_GPU_scGPU_50cut_5modes_scnoEV',
'benchmark_GPU_scorig_50cut_0modes',
'benchmark_GPU_scorig_50cut_1modes_2EV',
'benchmark_GPU_scorig_50cut_3modes',
'benchmark_GPU_scorig_50cut_5modes',
'benchmark_GPU_scorig_50cut_5modes_2EV',
'benchmark_GPU_scorig_50cut_5modes_3EV',
'benchmark_GPU_scorig_50cut_5modes_halfEV',
'benchmark_GPU_scorig_50cut_5modes_point7EV']
#BM5m_0_5ev = "benchmark_GPU_scorig_50cut_5modes_halfEV"

#a5BM.loadBenchmark( path_folder, BM5m )
#a5BM.loadBenchmark( path_folder, BM0m)
#a5BM.loadBenchmark( path_folder, BM3m)
#a5BM.loadBenchmark( path_folder, BM5_0_5mev)
#noModes = load_benchmarks( path_folder, name_benchmark_0modes )
#noModes_result = evaluate(noModes)


#n5Modes = load_benchmarks( path_folder, name_benchmark_5modes )
#n5Modes_result = evaluate(n5Modes)
#print n5Modes_result['1ACB'][idx_mean_10]

#min_10_ =  a5BM.getSorted('mean_10',BM5m, 0, True )
#min_sort_10 = a5BM.getSorted('mean_sort_10',BM5m,0, True )

#print "5 modes min 10 ", min_10, "5 modes min_10 sorted", min_sort_10


#nomod_min_10_ =  a5BM.getSorted('mean_10',BM0m, 0, True )
#nomod_min_sonrt_10 = a5BM.getSorted('mean_sort_10',BM0m,0, True )
#print "0 modes min 10 ", nomo_min_10, "0 modes min_10 sorted", nomod_min_sort_10

a5BM.plotRmsd(1000,BM0m,'1ACB' , use_energy= True)
a5BM.scoringPerformance( BM0m, '1ACB')
energies = a5BM.getDataProtein(BM5m, '1ACB', 'energy')

#plot scoring performance
plt.clf()
dict_scores = OrderedDict()
names = a5BM.getSorted('mean_10', BM5m,onlygetName = True)
for name in names:
    score = a5BM.scoringPerformance(BM5m, name)
    dict_scores[name] = score
xticks, y = zip(*dict_scores.items())
x = np.arange(len(xticks))

plt.xticks(x, xticks)
plt.ylim(-0.1, 1.1)
plt.plot(x, y, 'bo')

# plot rating of scored results
plt.clf()
dict_scores = OrderedDict()
names = a5BM.getSorted('mean_10', BM5m,onlygetName = True)
for name in names:
    score = a5BM.getWeightedResult(BM5m, name)
    if score is None:
        dict_scores[name] = -1
    else:
        dict_scores[name] = score
xticks,y = zip(*dict_scores.items())
x = np.arange(len(xticks))


plt.xticks( x, xticks)
plt.ylim(-0.1, 1.1)
plt.plot(x,y, 'bo')


for name in names:
    score = a5BM.getWeightedResult(BM0m, name)
    if score is None:
        dict_scores[name] = -1
    else:
        dict_scores[name] = score
xticks,y = zip(*dict_scores.items())
x = np.arange(len(xticks))


plt.xticks( x, xticks)
plt.ylim(-0.1, 1.1)
plt.plot(x,y, 'ro')


# plot rating of scored results TILLs method
plt.clf()
dict_scores = OrderedDict()
names = a5BM.getSorted('mean_10', BM5m,onlygetName = True)
for name in names:
    score = a5BM.getWeightedResultTill(BM5m, name)
    if score is None:
        dict_scores[name] = -1
    else:
        dict_scores[name] = score
xticks,y = zip(*dict_scores.items())
x = np.arange(len(xticks))


plt.xticks( x, xticks)
plt.ylim(-0.1, 1.1)
plt.plot(x,y, 'bo')


for name in names:
    score = a5BM.getWeightedResultTill(BM0m, name)
    if score is None:
        dict_scores[name] = -1
    else:
        dict_scores[name] = score
xticks,y = zip(*dict_scores.items())
x = np.arange(len(xticks))


plt.xticks( x, xticks)
plt.ylim(-0.1, 1.1)
plt.plot(x,y, 'ro')

#plot rmsd vs rmsd
plt.clf()
xData = np.zeros(len(names))
yData = np.zeros(len(names))
weight = 0
wcount = 0
for i,name in enumerate(names):
    x = a5BM.getDataProtein( BM0m, name, 'mean_10'  )
    y = a5BM.getDataProtein(BM5m, name, 'mean_10')
    xData[i] = x
    yData[i] = y
    if x is not None and x != -10000 and  y is not None and y != -10000 :
        wcount += 1
        weight += x/y
weight /= wcount

lin = np.arange(len(names))
plt.plot( lin, lin, 'r-')
plt.plot( xData,yData, 'bo')
plt.xlim(1,60)
plt.ylim(1,60)
print "ratio 0modes / 5modes", weight


#plot rmsd vs rmsd sorted
plt.clf()
xData = np.zeros(len(names))
yData = np.zeros(len(names))
weight = 0
wcount = 0
for i,name in enumerate(names):
    x = a5BM.getDataProtein( BM0m, name, 'mean_sort_10'  )
    y = a5BM.getDataProtein(BM5m, name, 'mean_sort_10')
    xData[i] = x
    yData[i] = y
    if x is not None and x != -10000 and  y is not None and y != -10000 :
        wcount += 1
        weight += x/y
weight /= wcount

lin = np.arange(len(names))
plt.plot( lin, lin, 'r-')
plt.plot( xData,yData, 'bo')
plt.xlim(1,60)
plt.ylim(1,60)
print "ratio 0modes / 5modes", weight, wcount

#plot rmsd vs rmsd sorted
plt.clf()
num_modes = 5
xData = np.zeros(len(names))
yData = np.zeros(len(names))
weight = 0
wcount = 0
for i,name in enumerate(names):
    x = a5BM.getDataProtein( BM5m, name, 'mean_10'  )
    y = a5BM.getDataProtein(BM5m, name, 'mean_sort_10')
    xData[i] = x
    yData[i] = y
    if x is not None and x != -10000 and  y is not None and y != -10000 :
        wcount += 1
        weight += x/y
weight /= wcount

lin = np.arange(len(names))
plt.plot( lin, lin, 'r-')
plt.plot( xData,yData, 'bo')
plt.xlim(1,60)
plt.ylim(1,60)
print "ratio mean rmsd /mean sorted rmsd", weight, wcount


#plot rmsd

a5BM.plotRmsd(30000, BM0m, '1ACB',use_energy= True)
print a5BM.getWeightedResultTill(BM0m, '1ACB')

#calc std_dev
#stdDevRmsd(self, name_benchmark, name_protein, maxIdx)

import pandas as pd

path_evaluation = "/home/glenn/Documents/Masterarbeit/analysis"
sort_BM = BM0m
BM = BM0m



for bm in BMLoad:
    df = pd.DataFrame(columns=['name', 'mean_10', 'rank_Till', 'rank', 'score', 'stdDev'])
    a5BM.loadBenchmark( path_folder, bm)
    names = a5BM.getSorted('mean_10', bm,onlygetName = True)
    row_list=[]
    for name in names:
        dict = {}
        dict['name'] = name
        dict['mean_10']    = a5BM.getDataProtein( bm, name, 'mean_10'  )
        dict['rank_Till'] = a5BM.getWeightedResultTill(bm, name)
        dict['rank']      = a5BM.getWeightedResult(bm, name)
        dict['score']      = a5BM.scoringPerformance(bm, name)
        dict['stdDev']      = a5BM.stdDevRmsd(bm, name, 50)
        #print dict
        row_list.append( dict )

    #df.append( dict, ignore_index=True )
    df = pd.DataFrame(row_list)
    file_name = os.path.join( path_evaluation, bm )
    df.to_csv(file_name, sep='\t', encoding='utf-8')

np.cov(np.asarray(frame['rank']), np.asarray( frame['rank_Till']))

correlation = np.cov(np.asarray(frame['rank']), np.asarray( frame['rank_Till']))/(np.std(np.asarray(frame['rank']), axis=0)*np.std(np.asarray(frame['rank_Till']), axis=0))
np.corrcoef(np.asarray(frame['rank']), np.asarray( frame['rank_Till']))