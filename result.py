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
path_folder= "/home/glenn/work/benchmark5_attract"

BMLoad = [
'benchmark_GPU_scorig_50cut_0modes',
'benchmark_GPU_scorig_nocut_1modes',
'benchmark_GPU_scorig_50cut_3modes',
'benchmark_GPU_scorig_50cut_5modes'
]

BMLoad = [
'dG_mr0_ml0_ev1p0_sO_c50_mr0_ml0_ev1',
'dG_mr10_ml10_ev1p0_sO_c50_mr10_ml10_ev1p0',
'dG_mr1_ml1_ev0p1_sO_c50_mr1_ml1_ev0p1',
'dG_mr1_ml1_ev0p5_sO_c50_mr1_ml1_ev0p5',
'dG_mr1_ml1_ev1p0_sO_c50_mr1_ml1_ev1p0',
'dG_mr1_ml1_ev2p0_sO_c50_mr1_ml1_ev2p0',
'dG_mr1_ml1_ev5p0_sO_c50_mr1_ml1_ev5p0',
'dG_mr3_ml3_ev1p0_sO_c50_mr3_ml3_ev1p0',
'dG_mr5_ml5_ev0p1_sO_c50_mr5_ml5_ev0p1',
'dG_mr5_ml5_ev0p5_sO_c50_mr5_ml5_ev0p5',
'dG_mr5_ml5_ev1p0_sO_c50_mr5_ml5_ev1',
'dG_mr5_ml5_ev2p0_sO_c50_mr5_ml5_ev2p0']
for bm in BMLoad:
    a5BM.loadBenchmark(path_folder, bm)


def evaluateModes1(path_folder, num_modes ):
    list_dir = [ d for d in os.listdir( path_folder ) if os.path.isdir(os.path.join( path_folder,d))]
    name_modefolder = "mode_eval"
    result = {}
    for name_prot in list_dir:
        dir_mode = os.path.join( os.path.join( path_folder, name_prot ),name_modefolder)
        frame_rec = pd.read_csv(os.path.join(dir_mode, "result_rec.dat"))
        frame_lig = pd.read_csv(os.path.join(dir_mode, "result_lig.dat"))
        rmsd_improve_rec = frame_rec['rmsd'][0] / frame_rec['rmsd'][num_modes]
        rmsd_improve_lig = frame_lig['rmsd'][0] / frame_lig['rmsd'][num_modes]
        mode_maxOverlap_rec = max(frame_rec['overlap'])
        mode_maxOverlap_lig = max(frame_lig['overlap'])
        result[name_prot] =[ rmsd_improve_rec, rmsd_improve_lig, mode_maxOverlap_rec, mode_maxOverlap_lig]
    return result

path_folder= "/home/glenn/work/benchmark5_attract"
mode_evaluation = evaluateModes1(path_folder, 19 )

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

data_dict ={}

path_evaluation = "/home/glenn/Documents/Masterarbeit/analysis/newrun"
#path_evaluation = "/home/glenn/Documents/Masterarbeit/analysis/Plots/0_5_modes_ORIGDOCK_ORIGSCORE"
deleted =['4GXU', '1DE4', '1EXB','3L89','4GAM','4FQI','1EER','1BGX','2HMI']
for bm in BMLoad:
    df = pd.DataFrame(columns=['name', 'mean_10', 'mean_sort_50','mean_sort_10','rank_Till', 'rank', 'score', 'stdDev', 'best_energy','num_10A','num_5A','num_sorted_10A','num_sorted_5A', 'num_10A_10', 'num_5A_10', 'num_10A_100', 'num_5A_100', 'num_10A_1000', 'num_5A_1000', 'one_star', 'two_star', 'three_star', 'num_g20A_20',"m_rmsd_impro_rec","m_rmsd_impro_lig", "m_overmax_rec", "m_overmax_lig"])
    #a5BM.loadBenchmark( path_folder, bm)
    names = a5BM.getSorted('mean_10', bm,onlygetName = True)
    names = sorted(names)

    row_list=[]
    for name in names:
        conti = False
        for del_name in deleted:
            if name == del_name:
                conti = True
        if conti:
            continue
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
            dict['num_g20A_20'] = a5BM.getNumGreater(bm, name, 20, maxindex=20, sorted=False)
            dict["m_rmsd_impro_rec"] = mode_evaluation[name][0]
            dict["m_rmsd_impro_lig"] =mode_evaluation[name][1]
            dict["m_overmax_rec"] =mode_evaluation[name][2]
            dict["m_overmax_lig"] =mode_evaluation[name][3]
        except:
            pass
        #print dict
        row_list.append( dict )

    #df.append( dict, ignore_index=True )
    df = pd.DataFrame(row_list)
    file_name = os.path.join( path_evaluation, bm )
    df.to_csv(file_name,  encoding='utf-8')
    data_dict[bm] =df

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

path_evaluation = "/home/glenn/Documents/Masterarbeit/analysis/newrun"
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




#9147/141/80143

#08912522171

######################GET NUM STARS OF ALL BENCHMARKS
list = []
for bm in benchmarks:
    dict = {}
    index_modes = bm.find("mr") + 2
    # index_modes_end = name_singleBench[0][index_modes:].find("_")
    num_modes = bm[index_modes:index_modes + 4].rsplit('_', 1)[0]
    dict['numModes']= int(num_modes)
    index_scale = bm.find("ev") + 2
    index_scale1 = bm.find("p") + 1
    scale = float(bm[index_scale]) + float(bm[index_scale1])/10
    dict['evScale']= scale
    dict['one_star_50'] =    benchmarks[bm]['one_star_50'].sum()
    dict['two_star_50'] = benchmarks[bm]['two_star_50'].sum()
    dict['three_star_50'] = benchmarks[bm]['three_star_50'].sum()
    list.append(dict)
frame = pd.DataFrame(list)
frame =  frame[['numModes','evScale','one_star_50','two_star_50','three_star_50']]
frame.to_csv("/home/glenn/test")


############################GET MAX DIFF OF STARS TO NO STARS
from collections import Counter
best_mode = 'dG_mr3_ml3_ev1p0_sO_c50_mr3_ml3_ev1p0'
list =[]
dict ={}
names = []
numpresent ={}
for bm in benchmarks:
    if bm !='dG_mr0_ml0_ev1p0_sO_c50_mr0_ml0_ev1':
        b = benchmarks[bm]
        diff_largest = (b['two_star_50']+b['three_star_50']- benchmarks[no_mode]['two_star_50']- benchmarks[no_mode]['three_star_50']).nlargest(10)

        names_largest = b.iloc[diff_largest.index]['name']
        diff_smallest = (b['two_star_50'] + b['three_star_50'] -  benchmarks[no_mode]['two_star_50'] -  benchmarks[no_mode][
            'three_star_50']).nsmallest(10)
        names_smallest = b.iloc[diff_smallest.index]['name']
        names = names + names_smallest.tolist()

        best_three  = b['three_star_50'].nlargest(3)
        best_three_name = b.iloc[best_three.index]
        best_two = b['two_star_50'].nlargest(3)
        best_two_name = b.iloc[best_two.index]

        print bm ,'\n', best_three_name[['name','three_star_50']], "\n",best_two_name[['name','two_star_50']],'\n','\n',diff_largest.values,'\n',names_largest.values, '\n','\n',diff_smallest.values,'\n',names_smallest.values,'\n\n\n'

counter =  Counter(names)
path ='/home/glenn/work/benchmark5_attract'

list = ['3MXW',
'1KXQ',
'3PC8',
'1IJK',
'1AK4',
'1RV6',
'1K74',
'1SYX',
'1K4C',
'4G6M',
'1XQS',
'1WEJ',
'1WDW',
'1QFW',
'2MTA',
'2OUL',
'1GL1',
'1BVK',
'2I25',
'1AY7',
'2YVJ',
'1FSK',
'1M27',
'1XD3',
'7CEI',
'2PCC',
'3VLB',
'1BUH',
'2CFH']
import shutil
for i in list:
    src = os.path.join(path,i)
    shutil.copytree(src, os.path.join('/home/glenn/work/benchmark5_best',i))


[
'1AK4A-refe.pdb',
'1K4CB-refe.pdb',
'1RV6A-refe.pdb',
'1XQSA-refe.pdb',
'4G6MA-refe.pdb',
'1AK4B-refe.pdb',
'1K74A-refe.pdb',
'1RV6B-refe2.pdb',
'3MXWA-refe.pdb',
'4G6MB-refe.pdb',
'1IJKA-refe.pdb',
'1K74B-refe.pdb',
'1RV6B-refe.pdb',
'3MXWB-refe.pdb',
'1IJKB-refe.pdb',
'1KXQA-refe.pdb',
'1SYXA-refe.pdb',
'3PC8A-refe.pdb',
'1K4CA-refe.pdb',
'1KXQB-refe.pdb',
'1SYXB-refe.pdb' ,
'3PC8B-refe.pdb']