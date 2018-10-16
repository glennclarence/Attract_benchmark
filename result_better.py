






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


path_folder= "/home/glenn/work/benchmark5_best"
path_evaluation="/home/glenn/Documents/Masterarbeit/analysis/181012_analysis"

path_csv =  "/home/glenn/Documents/Masterarbeit/analysis/181012_analysis"
path_base = "/home/glenn/Documents/Masterarbeit/analysis/181012_analysis/Plots"


BMLoad = []


file = open("/home/glenn/work/benchmarks")
for line in file.readlines():
    BMLoad.append(line.split()[0])




#
for bm in BMLoad:
    a5BM.loadBenchmark(path_folder, bm)
#



data_dict ={}


#path_evaluation = "/home/glenn/Documents/Masterarbeit/analysis/Plots/0_5_modes_ORIGDOCK_ORIGSCORE"
#deleted =['4GXU', '1DE4', '1EXB','3L89','4GAM','4FQI','1EER','1BGX','2HMI']
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
    bench = pd.read_csv(os.path.join( path_evaluation, bm))
    benchmarks[bm] = bench


list = []
for bm in benchmarks:
    dict = {}
    index_modesRec = bm.find("mr") + 2
    # index_modes_end = name_singleBench[0][index_modes:].find("_")
    num_modesRec = bm[index_modesRec:index_modesRec + 4].rsplit('_', 1)[0]
    dict['numModesRec']= int(num_modesRec)

    index_modesLig = bm.find("ml") + 2
    # index_modes_end = name_singleBench[0][index_modes:].find("_")
    num_modesLig = bm[index_modesLig:index_modesLig + 4].rsplit('_', 1)[0]
    dict['numModesLig'] = int(num_modesLig)


    index_scale = bm.find("ev") + 2
    index_scale1 = bm.find("p") + 1
    index_scale1
    scale = float(bm[index_scale:index_scale+2].replace("p","")) + float("0."+bm[index_scale1:index_scale1+3].replace("_","").replace("s",""))
    dict['evScale']= scale
    dict['one_star_50'] =    benchmarks[bm]['one_star_50'].sum()
    dict['two_star_50'] = benchmarks[bm]['two_star_50'].sum()
    dict['three_star_50'] = benchmarks[bm]['three_star_50'].sum()
    dict['sum_star'] = benchmarks[bm]['three_star_50'].sum()+benchmarks[bm]['one_star_50'].sum()+benchmarks[bm]['two_star_50'].sum()
    list.append(dict)
frame = pd.DataFrame(list)
frame =  frame[['numModesRec','numModesLig','evScale','one_star_50','two_star_50','three_star_50','sum_star']]
frame.to_csv( "/home/glenn/Documents/Masterarbeit/analysis/181012_analysis/scaleComparison/data.csv")




#names = a5BM.getSorted('mean_10', 'benchmark_ORIG_scorig_0modes',onlygetName = True)
names = []
#file = open('/home/glenn/work/benchmark5_attract/dir.txt')

#for name in file.readlines():
#    names.append(name.split()[0])
names=[]

dict_data = {}
for bm in BMLoad:
    frame = pd.read_csv(os.path.join(path_csv, bm))
    dict_data[bm] = frame


for name in names:
    path_protein = os.path.join(path_base, name)
    if not os.path.exists(path_protein):
        os.mkdir(path_protein)
    bm_mean = {}
    bm_star = {}
    bm_num = {}
    bm_rank = {}
    try:
        for bm in BMLoad:

            # plot =  a5BM.plotRmsd( max_indx = None, name_benchmark = bm, name_protein = name, use_energy = True, show = False, clear = False)
            readframe = pd.read_csv(os.path.join(path_base, name) + "/{}_rmsd.dat".format(bm))
            plt.clf()
            plt.scatter(readframe['rmsd'].values, readframe['energy '].values)
            plt.xlabel('RMSD (A)')
            plt.ylabel('Energy ')
            plt.title('{}\t{}'.format(bm, name))
            plt.ylim(ymax=1)
            plt.xlim(xmin=1, xmax=90)
            plt.savefig(os.path.join(path_protein, "fig-{}-{}.svg".format(name, bm)))
            # data = dict_data[bm]
            # dat = data.loc[data['name'] == name]
            # dict = {}
            # dict_mean = {}
            # dict_star = {}
            # dict_num = {}
            # dict_rank = {}
            # for column in dat:
            #     dict_mean['mean_10'] = dat['mean_10'].values[0]
            #     dict_mean['mean_sort_50'] = dat['mean_sort_50'].values[0]
            #     dict_mean['mean_sort_10'] = dat['mean_sort_10'].values[0]
            #     dict_mean['best_energy'] = dat['best_energy'].values[0]
            #     dict_star['one_star_50'] = dat['one_star_50'].values[0]
            #     dict_star['two_star_50'] = dat['two_star_50'].values[0]
            #     dict_star['three_star_50'] = dat['three_star_50'].values[0]
            #     dict_rank['rank_Till'] = dat['rank_Till'].values[0]
            #     dict_rank['rank'] = dat['rank'].values[0]
            #     dict_rank['score'] = dat['score'].values[0]
            #     dict_num['num_10A'] = dat['num_10A'].values[0]
            #     dict_num['num_5A'] = dat['num_5A'].values[0]
            #     dict_num['num_10A_10'] = dat['num_10A_10'].values[0]
            #     dict_num['num_5A_10'] = dat['num_5A_10'].values[0]
            #     dict_num['num_10A_100'] = dat['num_10A_100'].values[0]
            #     dict_num['num_5A_100'] = dat['num_5A_100'].values[0]
            #     dict_num['num_10A_1000'] = dat['num_10A_1000'].values[0]
            #     dict_num['num_5A_1000'] = dat['num_5A_1000'].values[0]
            #     dict_num['num_sorted_10A'] = dat['num_sorted_10A'].values[0]
            #     dict_num['num_sorted_5A'] = dat['num_sorted_5A'].values[0]
            # bm_mean[bm] = dict_mean
            # bm_star[bm] = dict_star
            # bm_num[bm] = dict_num
            # bm_rank[bm] = dict_rank

        # plt.clf()
        # for i, bm in enumerate(BMLoad):
        #     plt.bar(np.arange(len(bm_mean[bm])) * 2 * len(BMLoad) * 0.2 + i * 0.2, list(bm_mean[bm].values()), width=0.2,
        #             align='center')
        # plt.xticks(np.arange(len(bm_mean[bm])) * 2 * len(BMLoad) * 0.2, list(bm_mean[bm].keys()), rotation='vertical')
        # plt.legend()
        # plt.tight_layout()
        # plt.savefig(os.path.join(path_protein, "bar-mean-{}-{}.svg".format(name, bm)))
        #
        # plt.clf()
        # for i, bm in enumerate(BMLoad):
        #     plt.bar(np.arange(len(bm_rank[bm])) * 2 * len(BMLoad) * 0.2 + i * 0.2, list(bm_rank[bm].values()), width=0.2,
        #             align='center')
        # plt.xticks(np.arange(len(bm_rank[bm])) * 2 * len(BMLoad) * 0.2, list(bm_rank[bm].keys()), rotation='vertical')
        # plt.legend()
        # plt.tight_layout()
        # plt.savefig(os.path.join(path_protein, "bar-rank-{}-{}.svg".format(name, bm)))
        #
        # plt.clf()
        # for i, bm in enumerate(BMLoad):
        #     plt.bar(np.arange(len(bm_star[bm])) * 2 * len(BMLoad) * 0.2 + i * 0.2, list(bm_star[bm].values()), width=0.2,
        #             align='center')
        # plt.xticks(np.arange(len(bm_star[bm])) * 2 * len(BMLoad) * 0.2, list(bm_star[bm].keys()), rotation='vertical')
        # plt.legend()
        # plt.tight_layout()
        # plt.savefig(os.path.join(path_protein, "bar-star-{}-{}.svg".format(name, bm)))
        #
        # plt.clf()
        # for i, bm in enumerate(BMLoad):
        #     plt.bar(np.arange(len(bm_num[bm])) * 2 * len(BMLoad) * 0.2 + i * 0.2, list(bm_num[bm].values()), width=0.2,
        #             align='center')
        # plt.xticks(np.arange(len(bm_num[bm])) * 2 * len(BMLoad) * 0.2, list(bm_num[bm].keys()), rotation='vertical')
        # plt.legend()
        # plt.tight_layout()
        # plt.savefig(os.path.join(path_protein, "bar-num-{}-{}.svg".format(name, bm)))
    except:
        print "failed", name


for name in names:
    for bm in BMLoad:
        try:
            path_protein = os.path.join(path_base, name)
            if not os.path.exists(path_protein):
                os.mkdir(path_protein)
            file = open(os.path.join(path_base, name) + "/{}_rmsd.dat".format(bm), 'w+')
            file.write(",energy,rmsd \n")
            rmsd = a5BM._dict_dataBenchmark[bm][name][a5BM._dict_index['rmsd']][:]
            energy = a5BM._dict_dataBenchmark[bm][name][a5BM._dict_index['energy']][:]
            for i in range(len(a5BM._dict_dataBenchmark[bm][name][a5BM._dict_index['rmsd']])):
                file.write("{},{},{}\n".format(i, energy[i], rmsd[i]))
            file.close()
        except:
            print "filed", name
            pass

#plot rmsd vs rmsd
# plt.clf()
# xData = np.zeros(len(names))
# yData = np.zeros(len(names))
# weight = 0
# wcount = 0
# for i,name in enumerate(names):
#     x = a5BM.getDataProtein( BM0m, name, 'mean_10'  )
#     y = a5BM.getDataProtein(BM5m, name, 'mean_10')
#     xData[i] = x
#     yData[i] = y
#     if x is not None and x != -10000 and  y is not None and y != -10000 :
#         wcount += 1
#         weight += x/y
# weight /= wcount
#
# lin = np.arange(len(names))
# plt.plot( lin, lin, 'r-')
# plt.plot( xData,yData, 'bo')
# plt.xlim(1,60)
# plt.ylim(1,60)
# print "ratio 0modes / 5modes", weight
#
#
#
# dict_data = {}
# for bm in BMLoad:
#     frame = pd.read_csv(os.path.join(path_csv, bm))
#     dict_data[bm] = frame
#
#
# for i, bm in enumerate(BMLoad):
#     plt.bar(np.arange(len(bm_num[bm])) * 2 * len(BMLoad) * 0.2 + i * 0.2, list(bm_num[bm].values()), width=0.2,
#     align='center')
#     plt.xticks(np.arange(len(bm_num[bm])) * 2 * len(BMLoad) * 0.2, list(bm_num[bm].keys()), rotation='vertical')
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(os.path.join(path_protein, "bar-num-{}-{}.svg".format(name, bm)))
#
#
#

barbm = BMLoad
barnames = ['1AK4',
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

for i,name in enumerate(barnames):
    one = [];
    two = [];
    three = []
    for k,bm in enumerate(barbm):
        try:
            plt.gca().set_prop_cycle(None)
            one.append(dict_data[bm]['one_star_50'].loc[dict_data[bm]['name'] == name].values[0])
            two.append(dict_data[bm]['two_star_50'].loc[dict_data[bm]['name'] == name].values[0])
            three.append(dict_data[bm]['three_star_50'].loc[dict_data[bm]['name'] == name].values[0])
            #print name, bm, dict_data[bm]['one_star_50'].loc[dict_data[bm]['name'] == name].values[0],dict_data[bm]['two_star_50'].loc[dict_data[bm]['name'] == name].values[0],dict_data[bm]['three_star_50'].loc[dict_data[bm]['name'] == name].values[0]
            width = 1
            ind = np.arange(len(barbm))
            modes = [0, 3, 5, 10]
            p3 = plt.bar(ind +len(barbm)*i , one, width)
            p3 = plt.bar(ind+len(barbm)*i , two, width, bottom=one)
            p3 = plt.bar(ind+len(barbm)* i, three, width, bottom=two)
            plt.xticks(ind+len(barbm)* i, modes)
        except:
            pass


frame = pd.read_csv( "/home/glenn/Documents/Masterarbeit/analysis/181012_analysis/scaleComparison/data.csv")
numModesRec = [1,2,3,5,10,20]; numModesLig = [0,1]; maxScale = 10
sub = frame.loc[ (frame['numModesRec'].isin(numModesRec)) & (frame['numModesLig'].isin(numModesLig)) & (frame['evScale'] < maxScale)].sort_values(by=['evScale','numModesRec','numModesLig'])
#plt.scatter(sub['evScale'], sub['sum_star'])
#plt.show()
plt.clf()
one = sub['one_star_50'].tolist()
two = sub['two_star_50'].tolist()
three = sub['three_star_50'].tolist()
mrec = sub['numModesRec'].tolist()
mlig = sub['numModesLig'].tolist()
scale = sub['evScale'].tolist()
width = 1
ind = np.arange(len(one))
plt.bar(ind, one, width)
plt.bar(ind, two, width, bottom= one)
plt.bar(ind, three, width , bottom = two)
ticks = []
for a,b,c in zip(scale, mrec, mlig):
    ticks.append("{}\n{}\n{}".format(a,b,c))
plt.xticks(ind,ticks)
plt.show()
scale