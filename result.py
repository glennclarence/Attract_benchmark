import os
#from  load_pdbs import ProteinEsemble
import numpy as np
import matplotlib.pyplot as plt
from eval_class import ResultClass
from collections import OrderedDict
import pandas as pd

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
dict_indices['irmsd'] = 12
dict_indices['fnat'] = 13
dict_indices['num_10A_10'] = 14
dict_indices['num_5A_10'] =15
dict_indices['num_10A_100'] =16
dict_indices['num_5A_100'] =17
dict_indices['num_10A_1000']=18
dict_indices['num_5A_1000'] =19
dict_indices['one_star_50'] = 20
dict_indices['two_star_50'] = 21
dict_indices['three_star_50'] = 22
dict_indices['capri'] = 23


a5BM = ResultClass( dict_indices )

path_folder = "/home/glenn/Documents/Masterarbeit/test_bencheval/"
path_folder= "/home/glenn/cluster/benchmark5_attract"

BMLoad = [
'benchmark_GPU_scorig_50cut_0modes',
'benchmark_GPU_scorig_nocut_1modes',
'benchmark_GPU_scorig_50cut_3modes',
'benchmark_GPU_scorig_50cut_5modes'
]

for bm in BMLoad:
    a5BM.loadBenchmark(path_folder, bm)




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
    x = a5BM.getDataProtein( BM0m, name, 'mean_sort_10' )
    y = a5BM.getDataProtein( BM5m, name, 'mean_sort_10' )
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


sort_BM = BM0m
BM = BM0m


path_evaluation = "/home/glenn/Documents/Masterarbeit/analysis"
path_evaluation = "/home/glenn/Documents/Masterarbeit/analysis/Plots/0_5_modes_ORIGDOCK_ORIGSCORE"
for bm in BMLoad:
    df = pd.DataFrame(columns=['name', 'mean_10', 'mean_sort_50','mean_sort_10','rank_Till', 'rank', 'score', 'stdDev', 'best_energy','num_10A','num_5A','num_sorted_10A','num_sorted_5A', 'num_10A_10', 'num_5A_10', 'num_10A_100', 'num_5A_100', 'num_10A_1000', 'num_5A_1000', 'one_star', 'two_star', 'three_star'])
    #a5BM.loadBenchmark( path_folder, bm)
    names = a5BM.getSorted('mean_10', bm,onlygetName = True)
    row_list=[]
    for name in names:
        dict = {}
        try:
            dict['name'] = name
            dict['mean_10']    = a5BM.getDataProtein( bm, name, 'mean_10'  )

            dict['one_star_50'] = a5BM.getDataProtein(bm, name, 'one_star_50')
            dict['two_star_50'] = a5BM.getDataProtein(bm, name, 'two_star_50')
            dict['three_star_50'] = a5BM.getDataProtein(bm, name, 'three_star_50')

            dict['mean_sort_50'] = a5BM.getDataProtein(bm, name, 'mean_sort_50')
            dict['mean_sort_10'] = a5BM.getDataProtein(bm, name, 'mean_sort_10')
            dict['rank_Till'] = a5BM.getWeightedResultTill(bm, name)
            dict['rank']      = a5BM.getWeightedResult(bm, name)
            dict['score']      = a5BM.scoringPerformance(bm, name)
            dict['stdDev']      = a5BM.stdDevRmsd(bm, name, 50)
            dict['best_energy'] = min(a5BM.getDataProtein( bm, name, 'energy'  ))

            dict['num_10A'] = a5BM.getNumSmaller(bm, name, 10, sorted=False)
            dict['num_5A'] = a5BM.getNumSmaller(bm, name, 5, sorted=False)
            dict['num_10A_10'] = a5BM.getNumSmaller(bm, name, 10,maxindex=10,sorted = False)
            dict['num_5A_10'] = a5BM.getNumSmaller(bm, name, 5,maxindex=10,sorted = False)
            dict['num_10A_100'] = a5BM.getNumSmaller(bm, name, 10,maxindex=100, sorted=False)
            dict['num_5A_100'] = a5BM.getNumSmaller(bm, name, 5, maxindex=100,sorted=False)
            dict['num_10A_1000'] = a5BM.getNumSmaller(bm, name, 10,maxindex=1000, sorted=False)
            dict['num_5A_1000'] = a5BM.getNumSmaller(bm, name, 5,maxindex=1000, sorted=False)
            dict['num_sorted_10A'] = a5BM.getNumSmaller(bm, name, 10,sorted = True)
            dict['num_sorted_5A'] = a5BM.getNumSmaller(bm, name, 5,sorted = True)
        except:
            pass
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
for bm in bms:
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
bms=[
'1_benchmark_GPU_scGPU_50cut_0modes',
'1_benchmark_GPU_scGPU_50cut_3modes',
'1_benchmark_GPU_scGPU_50cut_3modes_sc2EV',
'1_benchmark_GPU_scGPU_50cut_3modes_scnoEV',
'1_benchmark_GPU_scGPU_50cut_5modes',
'1_benchmark_GPU_scGPU_50cut_5modes_sc2EV',
'1_benchmark_GPU_scGPU_50cut_5modes_scnoEV',
'1_benchmark_GPU_scorig_50cut_0modes',
'1_benchmark_GPU_scorig_50cut_1modes_2EV',
'1_benchmark_GPU_scorig_50cut_3modes',
'1_benchmark_GPU_scorig_50cut_5modes',
'1_benchmark_GPU_scorig_50cut_5modes_2EV',
'1_benchmark_GPU_scorig_50cut_5modes_3EV',
'1_benchmark_GPU_scorig_50cut_5modes_halfEV'  ]
#######################################################################plot rmsds according to ranking

bms = ['1_benchmark_GPU_scorig_50cut_0modes',
'1_benchmark_GPU_scorig_50cut_3modes',
'1_benchmark_GPU_scorig_50cut_5modes']
BM = [ "benchmark_GPU_scorig_50cut_0modes","benchmark_GPU_scorig_50cut_3modes","benchmark_GPU_scorig_50cut_5modes" ]
path = '/home/glenn/Documents/Masterarbeit/analysis/Plots/BM5_GPU_eval'


for bb in BM:
    a5BM.loadBenchmark(path_folder,bb)



for index in range(len(bms)):
    bm = bms[index]
    path1 = os.path.join( path, bm)
    df = pd.read_csv(path1, sep='\t')
    B = BM[index]
    column = 'rank'
    sorted = df.sort_values(by=column,ascending=False)[:20]

    for i,name in enumerate(sorted['name']):
        #plt.clf()
        #a5BM.plotRmsd(10000, B, name, use_energy=True)

        name_dir = "/home/glenn/{}_{}".format(B,column)


        if not os.path.exists(name_dir):
            os.makedirs(name_dir)
        # os.makedirs(name_dir, exist_ok=True)
        #plt.savefig('{}/{}_{}_{}.svg'.format(name_dir,i,name,column))
        filename = '{}/{}.dat'.format(name_dir,column)
        file = open(filename, 'w' )
       # sorted.reindex(np.arange(len(sorted.index)))
        #sorted.reset_index(np.arange(len(sorted.index)))
        sorted.reset_index(drop=True, inplace=True)
        print sorted
        for i, name in enumerate(sorted['name']):
            string = " {0:3}\t{1:6}\t{2:23.6g}\t{3:23.6g}\n".format(i, name, sorted[column][i], sorted['mean_10'][i] )
            file.write( string )



BMLoad = [
'benchmark_GPU_scorig_50cut_0modes',
'benchmark_GPU_scorig_nocut_1modes',
'benchmark_GPU_scorig_50cut_3modes',
'benchmark_GPU_scorig_50cut_5modes'
]


bm_dict={}
for bm in BMLoad:
    frame = pd.read_csv(os.path.join(path_evaluation, bm), sep='\t')
    dict ={}
    for column in frame:
        if column != 'Unnamed: 0':
            try:
                dict[column] = frame[column].mean()
            except:
                pass
    bm_dict[bm] = dict


###BARPLOT
plt.clf()
for i,bm in enumerate(BMLoad):
    plt.bar(np.arange(len(bm_dict[bm]))*2*len(BMLoad)*0.2+i*0.2, list(bm_dict[bm].values()),width=0.2, align='center')


plt.xticks(np.arange(len(bm_dict[bm]))*2*len(BMLoad)*0.2, list(bm_dict[bm].keys()), rotation='vertical')

    # librarie


plt.legend()

# Show graphic
plt.show()

########END BARPLOTT''''''########
path_base = "/home/glenn/Documents/Masterarbeit/analysis/Plots/evaluation"
path_csv = "/home/glenn/Documents/Masterarbeit/analysis"
names = a5BM.getSorted('mean_10', 'benchmark_GPU_scorig_50cut_0modes',onlygetName = True)

path_base = "/home/glenn/Documents/Masterarbeit/analysis/Plots/0_5_modes_ORIGDOCK_ORIGSCORE/evaluation"
path_csv = "/home/glenn/Documents/Masterarbeit/analysis/Plots/0_5_modes_ORIGDOCK_ORIGSCORE"
names = a5BM.getSorted('mean_10', 'benchmark_ORIG_scorig_0modes',onlygetName = True)
dict_data = {}
for bm in BMLoad:
    frame = pd.read_csv(os.path.join(path_csv, bm), sep='\t')
    dict_data[bm] = frame

for name in names:
    path_protein =os.path.join( path_base, name)
    if not os.path.exists(path_protein):
        os.mkdir(path_protein)
    bm_mean = {}
    bm_star = {}
    bm_num = {}
    bm_rank = {}
    for bm in BMLoad:
        plot =  a5BM.plotRmsd( max_indx = None, name_benchmark = bm, name_protein = name, use_energy = True, show = False, clear = True)
        plot.savefig(os.path.join(path_protein, "fig-{}-{}.svg".format(name,bm)))
        data = dict_data[bm]
        dat =  data.loc[data['name'] == name]
        dict = {}
        dict_mean ={}
        dict_star = {}
        dict_num= {}
        dict_rank = {}
        for column in dat:
            dict_mean['mean_10'] = dat['mean_10'].values
            dict_mean['mean_sort_50'] = dat['mean_sort_50'].values
            dict_mean['mean_sort_10'] = dat['mean_sort_10'].values
            dict_mean['best_energy'] = dat['best_energy'].values
            dict_star['one_star_50'] = dat['one_star_50'].values
            dict_star['two_star_50'] = dat['two_star_50'].values
            dict_star['three_star_50'] = dat['three_star_50'].values
            dict_rank['rank_Till'] = dat['rank_Till'].values
            dict_rank['rank'] = dat['rank'].values
            dict_rank['score'] = dat['score'].values
            dict_num['num_10A'] = dat['num_10A'].values
            dict_num['num_5A'] = dat['num_5A'].values
            dict_num['num_10A_10'] = dat['num_10A_10'].values
            dict_num['num_5A_10'] = dat['num_5A_10'].values
            dict_num['num_10A_100'] = dat['num_10A_100'].values
            dict_num['num_5A_100'] = dat['num_5A_100'].values
            dict_num['num_10A_1000'] = dat['num_10A_1000'].values
            dict_num['num_5A_1000'] = dat['num_5A_1000'].values
            dict_num['num_sorted_10A'] = dat['num_sorted_10A'].values
            dict_num['num_sorted_5A'] = dat['num_sorted_5A'].values
        bm_mean[bm] = dict_mean
        bm_star[bm] = dict_star
        bm_num[bm] = dict_num
        bm_rank[bm] = dict_rank
    plt.clf()
    for i, bm in enumerate(BMLoad):
        plt.bar(np.arange(len(bm_mean[bm])) * 2 * len(BMLoad) * 0.2 + i * 0.2, list(bm_mean[bm].values()), width=0.2,
                align='center')
    plt.xticks(np.arange(len(bm_mean[bm])) * 2 * len(BMLoad) * 0.2, list(bm_mean[bm].keys()), rotation='vertical')      
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(path_protein, "bar-mean-{}-{}.svg".format(name,bm)))

    plt.clf()
    for i, bm in enumerate(BMLoad):
        plt.bar(np.arange(len(bm_rank[bm])) * 2 * len(BMLoad) * 0.2 + i * 0.2, list(bm_rank[bm].values()), width=0.2,
                align='center')
    plt.xticks(np.arange(len(bm_rank[bm])) * 2 * len(BMLoad) * 0.2, list(bm_rank[bm].keys()), rotation='vertical')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(path_protein, "bar-rank-{}-{}.svg".format(name, bm)))

    plt.clf()
    for i, bm in enumerate(BMLoad):
        plt.bar(np.arange(len(bm_star[bm])) * 2 * len(BMLoad) * 0.2 + i * 0.2, list(bm_star[bm].values()), width=0.2,
                align='center')
    plt.xticks(np.arange(len(bm_star[bm])) * 2 * len(BMLoad) * 0.2, list(bm_star[bm].keys()), rotation='vertical')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(path_protein, "bar-star-{}-{}.svg".format(name, bm)))

    plt.clf()
    for i, bm in enumerate(BMLoad):
        plt.bar(np.arange(len(bm_num[bm])) * 2 * len(BMLoad) * 0.2 + i * 0.2, list(bm_num[bm].values()), width=0.2,
                align='center')
    plt.xticks(np.arange(len(bm_num[bm])) * 2 * len(BMLoad) * 0.2, list(bm_num[bm].keys()), rotation='vertical')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(path_protein, "bar-num-{}-{}.svg".format(name, bm)))




