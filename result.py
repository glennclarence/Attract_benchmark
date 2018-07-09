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
    df = pd.DataFrame(columns=['name', 'mean_10', 'mean_sort_50','rank_Till', 'rank', 'score', 'stdDev', 'best_energy'])
    a5BM.loadBenchmark( path_folder, bm)
    names = a5BM.getSorted('mean_10', bm,onlygetName = True)
    row_list=[]
    for name in names:
        dict = {}
        dict['name'] = name
        dict['mean_10']    = a5BM.getDataProtein( bm, name, 'mean_10'  )
        dict['mean_sort_50'] = a5BM.getDataProtein(bm, name, 'mean_sort_50')
        dict['rank_Till'] = a5BM.getWeightedResultTill(bm, name)
        dict['rank']      = a5BM.getWeightedResult(bm, name)
        dict['score']      = a5BM.scoringPerformance(bm, name)
        dict['stdDev']      = a5BM.stdDevRmsd(bm, name, 50)
        dict['best_energy'] = min(a5BM.getDataProtein( bm, name, 'energy'  ))
        #print dict
        row_list.append( dict )

    #df.append( dict, ignore_index=True )
    df = pd.DataFrame(row_list)
    file_name = os.path.join( path_evaluation, bm )
    df.to_csv(file_name, sep='\t', encoding='utf-8')

np.cov(np.asarray(frame['rank']), np.asarray( frame['rank_Till']))

correlation = np.cov(np.asarray(frame['rank']), np.asarray( frame['rank_Till']))/(np.std(np.asarray(frame['rank']), axis=0)*np.std(np.asarray(frame['rank_Till']), axis=0))
np.corrcoef(np.asarray(frame['rank']), np.asarray( frame['rank_Till']))




BMLoad = [
'1_benchmark_GPU_scorig_50cut_0modes',
'1_benchmark_GPU_scorig_50cut_1modes_2EV',
'1_benchmark_GPU_scorig_50cut_3modes',
'1_benchmark_GPU_scorig_50cut_5modes',
'1_benchmark_GPU_scorig_50cut_5modes_2EV',
'1_benchmark_GPU_scorig_50cut_5modes_3EV',
'1_benchmark_GPU_scorig_50cut_5modes_halfEV']

path_evaluation = "/home/glenn/Documents/Masterarbeit/analysis"
benchmarks = {}
for bm in BMLoad:
    bench = pd.read_csv(os.path.join( path_evaluation, bm), sep = "\t")
    benchmarks[bm] = bench

plt.clf()
num_plots = 14

# Have a look at the colormaps here and decide which one you'd like:
# http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
colormap = plt.cm.gist_ncar
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])
for name, data in benchmarks.iteritems():
    plt.plot( data['mean_10'], 'o',label=name)
plt.legend(framealpha=0.5)

name_Plotpath = "/home/glenn/Documents/Masterarbeit/analysis/Plots"
name_plot = "1_benchmark_GPU_scorig_50cut_5modes_halfEV"
import pandas as pd


def loadPlotData(name_path, name_plot):
    data = pd.read_csv(os.pat.join(name_path, name_plot), sep="\t")
    return data


def readTiming( name_path, name_plot):
    timing = pd.read_csv(os.path.join( name_path, name_plot), sep="\s+")
    return timing

t5mgpu = readTiming( "/home/glenn/cluster", "timing_GPU_5modes")
t0mgpu = readTiming( "/home/glenn/cluster", "timing_GPU_0modes")

dict0ml = []
for protein in t0morig['Protein']:
    print protein
    bla = t0mgpu.loc[t0mgpu['Protein'] == protein, 'time'].tolist()
    dict0ml.append( (protein, bla[0]) )

dict5ml = []
for protein in t5morig['Protein']:
    print protein
    bla = t5mgpu.loc[t5mgpu['Protein'] == protein, 'time'].tolist()
    dict5ml.append( (protein, float(bla[0])) )
df5mgpu = pd.DataFrame(dict5ml, dtype = 'float',columns = ['Protein', 'time'])
speedup = np.asarray(t5morig['time']/df5mgpu['time'] )

path = '/home/glenn/cluster/benchmark5_attract/'
this_list = []
for protein in df5mgpu['Protein']:
    filename_dofs=path + protein +'/input/dofs.dat'
    #print filename_dofs
    num = 0
    idx= 1
    with open(filename_dofs, 'r') as f:
        for line in f.readlines():
            if line.startswith("#{}".format(idx)) or line.startswith("# {}".format(idx)):
                #print line
                num = idx
                idx += 1

    filename_pdb = path + protein + '/receptor-for-docking.pdb'
    filename_pdb = path + protein + '/input/' + protein + "-receptor-for-docking-reduce.pdb"
    with open(filename_pdb, 'r') as f:
       len_rec = len( f.readlines())
    filename_pdb = path + protein + '/ligand-for-docking.pdb'
    filename_pdb = path + protein + '/input/'+protein+"-ligand-for-docking-reduce.pdb"
    with open(filename_pdb, 'r') as f:
        len_lig = len(f.readlines())
    this_list.append((protein, num, len_rec, len_lig))


this_list = []
t5mgpu = readTiming( "/home/glenn/cluster", "timing_GPU_5modes")
for protein in t5mgpu['Protein']:
    filename_dofs=path + protein +'/input/dofs.dat'
    #print filename_dofs
    num = 0
    idx= 1
    with open(filename_dofs, 'r') as f:
        for line in f.readlines():
            if line.startswith("#{}".format(idx)) or line.startswith("# {}".format(idx)):
                #print line
                num = idx
                idx += 1

    filename_pdb = path + protein + '/receptor-for-docking.pdb'
    filename_pdb = path + protein + '/input/' + protein + "-receptor-for-docking-reduce.pdb"
    with open(filename_pdb, 'r') as f:
       len_rec = len( f.readlines())
    filename_pdb = path + protein + '/ligand-for-docking.pdb'
    filename_pdb = path + protein + '/input/'+protein+"-ligand-for-docking-reduce.pdb"
    with open(filename_pdb, 'r') as f:
        len_lig = len(f.readlines())
    this_list.append((protein, num, len_rec, len_lig))


1BGX_times=[
(5000,   32.0540018082),
(10000,  55.67700),
(20000,  105.064357),
(30000,  152.293265),
(40000,  203.344777),
(50000,  259.515478),
(60000,  302.787827),
(70000,  353.170102),
(80000,  402.038713)
]
x=[5000,
10000,
20000,
30000,
40000,
50000,
60000,
70000,
80000
]
y=[
32.0540018082,
55.67700,
105.064357,
152.293265,
203.344777,
259.515478,
302.787827,
353.170102,
402.03871]

np.polyfit(x,y,1)
array([  4.95791618e-03,   6.25623464e+00])
