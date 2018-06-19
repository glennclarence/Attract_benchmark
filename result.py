import os
from  load_pdbs import ProteinEsemble
import numpy as np
import matplotlib.pyplot as plt








path_folder = "/home/glenn/Documents/Masterarbeit/test_bencheval"

def load_benchmarks( path_folder):
    bench = list()
    for name_prot in os.listdir( path_folder ):
        dir_prot = os.path.join( path_folder, name_prot )
        benchmarks = [( d, os.path.join( dir_prot, d )) for d in os.listdir(dir_prot) if os.path.isdir( os.path.join( dir_prot,d)) and d != 'input' ]
        input =[ os.path.join(dir_prot, d) for d in os.listdir(dir_prot) if os.path.isdir(os.path.join(dir_prot, d)) and d == 'input']
        bench.append( (name_prot,input,benchmarks))
    return bench
ext_rmsd = "-rmsd.result"
ext_scored = "-sorted-dr.dat"
filename_pdb = "-receptor-for-docking"

def mode_ext( num_modes):
    return "-" + num_modes +"-modes.dat"

bechmarks = load_benchmarks( path_folder )



result = list()
for name_protein, input, name_bench in bechmarks:
    print name_protein
    result_protein = list()
    for name_singleBench in name_bench:
        index_modes = 27
        num_modes = name_singleBench[0][index_modes]
        print "\t",name_singleBench[0]
        if name_singleBench[0] != 'input':
            filename_rmsd = os.path.join(name_singleBench[1],name_protein+filename_pdb+ext_rmsd)
            with open(filename_rmsd, 'r' ) as f:
                lines = f.readlines()
                rmsd = np.zeros( len(lines), dtype = float)
                pos = np.arange( len(lines), dtype = int  )
                for i, line in enumerate(lines):
                    rmsd[i] = float( line.split()[-1])


            plt.plot( rmsd[:20], pos[:20])
            #plt.show()
            print "\t\t", rmsd.mean(), rmsd[:10].mean(), rmsd[:50].mean(), np.sort(rmsd)[:10].mean(), np.sort(rmsd)[:50].mean()


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
                    print "\t\t", "amplitudes", amplitude
                    print "\t\t", "eigenvalues", eigenvalues
        result_protein.append( (int(num_modes), rmsd, pos,  min(rmsd), rmsd.mean(), rmsd[:10].mean(), rmsd[:50].mean(), np.sort(rmsd)[:10].mean(), np.sort(rmsd)[:50].mean(),  amplitude, eigenvalues))
    result.append( (name_protein, result_protein) )


print result


result_rmsd= np.zeros( 10, dtype=float)
for protein , protein_result in result:
    for benchmark in protein_result:
        result_rmsd[benchmark[0]] += benchmark[5]

for i in range(len(result_rmsd)):
    result_rmsd[i] /= len(result)
    print "modes " , i, " ",result_rmsd[i]