




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
#'dG_mr1_ml1_ev1p0_sO_c50_mr1_ml1_ev1p0',
'dG_mr1_ml1_ev2p0_sO_c50_mr1_ml1_ev2p0',
'dG_mr1_ml1_ev5p0_sO_c50_mr1_ml1_ev5p0',
'dG_mr3_ml3_ev1p0_sO_c50_mr3_ml3_ev1p0',
'dG_mr5_ml5_ev0p1_sO_c50_mr5_ml5_ev0p1',
'dG_mr5_ml5_ev0p5_sO_c50_mr5_ml5_ev0p5',
'dG_mr5_ml5_ev1p0_sO_c50_mr5_ml5_ev1',
'dG_mr5_ml5_ev2p0_sO_c50_mr5_ml5_ev2p0']


for bm in BMLoad:
    a5BM.loadBenchmark(path_folder, bm)




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




benchmarks = {}
for bm in BMLoad:
    bench = pd.read_csv(os.path.join( path_evaluation, bm), sep = "\t")
    benchmarks[bm] = bench


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
    dict['ecScale']= scale
    dict['one_star_50'] =    benchmarks[bm]['one_star_50'].sum()
    dict['two_star_50'] = benchmarks[bm]['two_star_50'].sum()
    dict['three_star_50'] = benchmarks[bm]['three_star_50'].sum()
    list.append(dict)
frame = pd.DataFrame(list)
frame.to_csv("/home/glenn/test")