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
        bench.append( (name_prot,benchmarks))
    return bench
ext_rmsd = "-rmsd.result"
ext_scored = "-sorted-dr.dat"
filename_pdb = "-receptor-for-docking"
def mode_ext( num_modes):
    return "-" + num_modes +"-modes.dat"

bechmarks = load_benchmarks( path_folder )

for name_protein, name_bench in bechmarks:

    for name_singleBench in name_bench:
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
            print name_protein,name_singleBench[0],rmsd.mean(), rmsd[:10].mean(), rmsd[:50].mean(), np.sort(rmsd)[:10].mean(), np.sort(rmsd)[:50].mean()

            index_modes = 26
            num_modes = name_singleBench[0][index_modes]
           # filename_modes  = os.path.join(name_singleBench[1],name_protein+filename_pdb+modes_ext(num_modes))
            print filename_modes

