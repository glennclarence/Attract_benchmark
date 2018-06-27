import os
#from  load_pdbs import ProteinEsemble
import numpy as np
import matplotlib.pyplot as plt
from eval_class import ResultClass


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

attract5BM = ResultClass( dict_indices )

path_folder = "/home/glenn/Documents/Masterarbeit/test_bencheval/"
#path_folder= "/home/glenn/cluster/benchmark5_attract"

name_benchmark_0modes = "benchmark_GPU_scorig_50cut_0modes"

attract5BM.loadBenchmark( path_folder, name_benchmark_0modes)
attract5BM.loadBenchmark( path_folder, name_benchmark_0modes)
#noModes = load_benchmarks( path_folder, name_benchmark_0modes )
#noModes_result = evaluate(noModes)

name_benchmark_5modes = "benchmark_GPU_scorig_50cut_5modes"
#n5Modes = load_benchmarks( path_folder, name_benchmark_5modes )
#n5Modes_result = evaluate(n5Modes)
#print n5Modes_result['1ACB'][idx_mean_10]

print attract5BM.getMaximum('rmsd','benchmark_GPU_scorig_50cut_0modes',0, True )
print attract5BM.getSorted('rmsd','benchmark_GPU_scorig_50cut_0modes',0, True )


attract5BM.getDataProtein('benchmark_GPU_scorig_50cut_0modes',attract5BM.getMaximum('rmsd','benchmark_GPU_scorig_50cut_0modes',0 ,True))
