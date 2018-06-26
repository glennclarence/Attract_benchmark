import os
#from  load_pdbs import ProteinEsemble
import numpy as np
import matplotlib.pyplot as plt


idx_mode = 0
idx_rmsd = 1
idx_energy = 2
idx_pos = 3
idx_min = 4
idx_mean = 5
idx_mean_10 = 6
idx_mean_50 = 7
idx_mean_sort_10 = 8
idx_mean_sort_50 = 9
idx_amp = 10
idx_ev = 11
tot_num = 12

dict_indices = {}
dict_indices['idx_mode'] = 0
dict_indices['idx_rmsd'] = 1
dict_indices['idx_energy'] = 2
dict_indices['idx_pos'] = 3
dict_indices['idx_min'] = 4
dict_indices['idx_mean'] = 5
dict_indices['idx_mean_10'] = 6
dict_indices['idx_mean_50'] = 7
dict_indices['idx_mean_sort_10'] = 8
dict_indices['idx_mean_sort_50'] = 9
dict_indices['idx_amp'] = 10
dict_indices['idx_ev'] = 11

def getEnergyfromFile( filename_scoring):
    energies = {}
    with open(filename_scoring, 'r') as f:
        count = 0
        idx = 1
        for line in f.readlines():
            if line.startswith("#{}".format(idx)):
                idx += 1
                number =  float(line[1])
                count += 1
            elif line.startswith("## Energy:"):
               # print line
                energy = float(line.split()[2])
                count += 1
            if count == 2:
                energies[number] = energy
                count = 0
    return energies





def load_benchmarks( path_folder, name_benchmark):
    bench = list()
    list_dir = [ d for d in os.listdir( path_folder ) if os.path.isdir(os.path.join( path_folder,d))]
    for name_prot in list_dir:
        dir_prot = os.path.join( path_folder, name_prot )
        benchmarks = [( d, os.path.join( dir_prot, d )) for d in os.listdir(dir_prot) if os.path.isdir( os.path.join( dir_prot,d)) and d != 'input' and d == name_benchmark]
        input =[ os.path.join(dir_prot, d) for d in os.listdir(dir_prot) if os.path.isdir(os.path.join(dir_prot, d)) and d == 'input']
        bench.append( (name_prot,input,benchmarks))
    return bench
ext_rmsd = "-rmsd.result"
ext_scored = "-sorted-dr.dat"
filename_pdb = "-receptor-for-docking"

def mode_ext( num_modes):
    return "-" + num_modes +"-modes.dat"



def evaluate( bechmarks):
    result = {}
    count = 0
    for name_protein, input, name_bench in bechmarks:
        count += 1
        print count,name_bench[0][0], name_protein
        result_protein = [None]*tot_num
        for name_singleBench in name_bench:
            index_modes = name_singleBench[0].find("modes")-1
            num_modes = name_singleBench[0][index_modes]
            #print "\t",name_singleBench[0]
            if name_singleBench[0] != 'input':
                filename_rmsd = os.path.join(name_singleBench[1],name_protein+filename_pdb+ext_rmsd)
                with open(filename_rmsd, 'r' ) as f:
                    lines = f.readlines()
                    len_list = len(lines)
                    if len_list ==0:
                        continue
                    rmsd = np.zeros( len(lines), dtype = float)
                    pos = np.arange( len(lines), dtype = int  )
                    for i, line in enumerate(lines):
                        rmsd[i] = float( line.split()[-1])


                plt.plot( rmsd[:20], pos[:20])
                #plt.show()
                #print "\t\t", "total rmsd", rmsd.mean(), "\n\t\t10 rmsd",rmsd[:10].mean(), "\n\t\t50 rmsd",rmsd[:50].mean(), "\n\t\t10 sorted rmsd",np.sort(rmsd)[:10].mean(),"\n\t\t50 sorted rmsd", np.sort(rmsd)[:50].mean()


                if int(num_modes) > 0:
                    filename_modes  = os.path.join(input[0],name_protein+filename_pdb+mode_ext(num_modes))
                    with open(filename_modes, 'r') as f:
                        lines = f.readlines()
                        mode_idx = 1
                        eigenvalues = np.zeros( int(num_modes), dtype = float)
                        amplitude = np.zeros( int(num_modes) ,  dtype = float)
                        for i, line in enumerate(lines):
                            if len( line.split()) == 2 and int(line.split()[0]) == mode_idx:
                                mode_idx += 1
                                eigenvalues[mode_idx-2] = float(line.split()[1])
                            elif len(line.split()) == 3:
                                amplitude[mode_idx-2] += float(line.split()[0])*float(line.split()[0]) + float(line.split()[1])*float(line.split()[1]) + float(line.split()[2])*float(line.split()[2])
                        #print "\t\t", "amplitudes", amplitude
                        #print "\t\t", "eigenvalues", eigenvalues
                else:
                    amplitude = 0
                    eigenvalues = 0
                filename_scoring=os.path.join(name_singleBench[1],name_protein+filename_pdb+ext_scored)
                energies = getEnergyfromFile(filename_scoring)
                #print "energies",energies
            result_protein[idx_mode]    = int(num_modes)
            result_protein[idx_rmsd] = rmsd
            result_protein[idx_energy] = energies
            result_protein[idx_pos] = pos
            result_protein[idx_min] = min(rmsd)
            result_protein[idx_mean] = rmsd.mean()
            result_protein[idx_mean_10] = rmsd[:10].mean()
            result_protein[idx_mean_50] = rmsd[:50].mean()
            result_protein[idx_mean_sort_10] = np.sort(rmsd)[:10].mean()
            result_protein[idx_mean_sort_50] = np.sort(rmsd)[:50].mean()
            result_protein[idx_amp] = amplitude
            result_protein[idx_ev] = eigenvalues
           # result_protein.append( (int(num_modes), rmsd, energies,pos,  min(rmsd), rmsd.mean(), rmsd[:10].mean(), rmsd[:50].mean(), np.sort(rmsd)[:10].mean(), np.sort(rmsd)[:50].mean(),  amplitude, eigenvalues ))


        result[name_protein] = result_protein
    return result




#result_rmsd= np.zeros( 10, dtype=float)
#for protein , protein_result in result:
#    for benchmark in protein_result:
 #       result_rmsd[benchmark[0]] += benchmark[5]

#for i in range(len(result_rmsd)):
#    result_rmsd[i] /= len(result)
   # print "modes " , i, " ",result_rmsd[i]


#filename_scoring = "/home/glenn/cluster/benchmark_attract_test/1AVX/benchmark_ORI_scorig_50cut_5modes_2/1AVX-receptor-for-docking-sorted-dr.dat"
#energies =  getEnergyfromFile( filename_scoring)
#print energies


path_folder = "/home/glenn/Documents/Masterarbeit/test_bencheval/"
path_folder= "/home/glenn/cluster/benchmark5_attract"

name_benchmark_0modes = "benchmark_GPU_scorig_50cut_0modes"
noModes = load_benchmarks( path_folder, name_benchmark_0modes )
noModes_result = evaluate(noModes)

name_benchmark_5modes = "benchmark_GPU_scorig_50cut_5modes"
n5Modes = load_benchmarks( path_folder, name_benchmark_5modes )
n5Modes_result = evaluate(n5Modes)
print n5Modes_result['1ACB'][idx_mean_10]