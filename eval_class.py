
import os
#from  load_pdbs import ProteinEsemble
import numpy as np
import matplotlib.pyplot as plt

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



def evaluate( bechmarks, dict_indices):
    result = {}
    count = 0
    for name_protein, input, name_bench in bechmarks:
        count += 1
        print count,name_bench[0][0], name_protein
        result_protein = [None]*len(dict_indices)
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


                #plt.plot( rmsd[:20], pos[:20])
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
            result_protein[dict_indices['mode']]    = int(num_modes)
            result_protein[dict_indices['rmsd']] = rmsd
            result_protein[dict_indices['energy']] = energies
            result_protein[dict_indices['pos']] = pos
            result_protein[dict_indices['min']] = min(rmsd)
            result_protein[dict_indices['mean']] = rmsd.mean()
            result_protein[dict_indices['mean_10']] = rmsd[:10].mean()
            result_protein[dict_indices['mean_50']] = rmsd[:50].mean()
            result_protein[dict_indices['mean_sort_10']] = np.sort(rmsd)[:10].mean()
            result_protein[dict_indices['mean_sort_50']] = np.sort(rmsd)[:50].mean()
            result_protein[dict_indices['amp']] = amplitude
            result_protein[dict_indices['ev']] = eigenvalues
           # result_protein.append( (int(num_modes), rmsd, energies,pos,  min(rmsd), rmsd.mean(), rmsd[:10].mean(), rmsd[:50].mean(), np.sort(rmsd)[:10].mean(), np.sort(rmsd)[:50].mean(),  amplitude, eigenvalues ))


        result[name_protein] = result_protein
    return result


class ResultClass:
    def __init__(self, dict_index):

        self._dict_index = dict_index
        self._dict_dataBenchmark = {}
        self._activeBM = None
        self._activeBMReference = None
        self._activeProtein = None
        self._proteinNames = set()

    def add_dataBenchmark(self, name_benchmark, data_benchmark):
        self._dict_dataBenchmark[name_benchmark] = data_benchmark
        for name_protein, data in data_benchmark.iteritems():
            self._proteinNames.add(name_protein)

    def getNamesBenchmarks(self):
        return [name for name, value in self._dict_dataBenchmark.iteritems()]

    def getSortBenchmark(self, name_benchmark, idx_name):
        bm = self._dict_dataBenchmark[name_benchmark]
        return sorted(bm[idx_name],key=bm.get)

    def getDataProtein(self, name_benchmark, name_protein):
        return self._dict_dataBenchmark[name_benchmark][name_protein]

    def getDataProtein(self, name_benchmark, name_protein, name_index):
        return self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index[name_index]]

    def setActiveBenchmarks(self, name_benchmark, name_benchmarkReference =None):
        self._activeBM = name_benchmark
        self._activeBMReference = name_benchmarkReference

    def setActiveProtein(self, name_benchmark):
        self._activeBM = name_benchmark

    def setActiveProtein(self, name_protein):
        self._activeProtein = name_protein

    def plotRmsd(self, name_benchmark = None, name_protein = None, use_energy = False):
        if name_protein is None:
            name_protein = self._activeProtein
        if name_benchmark is None:
            name_benchmark = self._activeBM
        rmsd = self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index['rmsd']]
        if not use_energy:
            y = self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index['pos']]
        else:
            y = self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index['energy']]
        plt.plot( rmsd, y)
        plt.show()

    def getProteinNames(self):
        return self._proteinNames

    def getMaximum(self, name_index, name_benchmark = None, index_subarray = None, onlygetName = False):
        if name_benchmark is None:
            name_benchmark = self._activeBM
        dict_data = self._dict_dataBenchmark[name_benchmark]
        if not onlygetName:
            if index_subarray is None:
                return  min( dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]]))
            else:
                return  min( dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]][index_subarray]))
        if onlygetName:
            if index_subarray is None:
                return  max( dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]]))[0]
            else:
                return  max(dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]][index_subarray]))[0]

    def getMaximum(self, name_index, name_benchmark = None, index_subarray = None, onlygetName = False):
        if name_benchmark is None:
            name_benchmark = self._activeBM
        dict_data = self._dict_dataBenchmark[name_benchmark]
        if not  onlygetName:
            if index_subarray is None:
                return  max( dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]]))
            else:
                return  max(dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]][index_subarray]))

        if onlygetName:
            if index_subarray is None:
                return  max( dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]]))[0]
            else:
                return  max(dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]][index_subarray]))[0]

    def getSorted(self, name_index, name_benchmark=None, index_subarray=None, onlygetName = False):
        if name_benchmark is None:
            name_benchmark = self._activeBM
        dict_data = self._dict_dataBenchmark[name_benchmark]
        if index_subarray is None:
            if onlygetName:
                return [d[0] for d in sorted(dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]]))]
            else:
                return sorted(dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]]))
        else:
            if onlygetName:
                return [d[0] for d in sorted(dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]][index_subarray]))]
            else:
                return  sorted(dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]][index_subarray]))

    def loadBenchmark(self , path_benchmark, name_benchmark):
        data_benchmark = load_benchmarks(path_benchmark, name_benchmark)
        data_evaluated = evaluate(data_benchmark, self._dict_index)
        self.add_dataBenchmark( name_benchmark, data_evaluated)
