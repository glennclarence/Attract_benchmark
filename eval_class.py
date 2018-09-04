
import os
#from  load_pdbs import ProteinEsemble
import numpy as np
import matplotlib.pyplot as plt

def getEnergyfromFile( filename_scoring):
    energies = list()
    #try:
    with open(filename_scoring, 'r') as f:
        count = 0
        idx = 1
        for line in f.readlines():
            if line.startswith("#{}".format(idx)):

                number =  float(line[1])
                idx += 1
                count += 1
            elif line.startswith("## Energy:"):
               # print line
                energy = float(line.split()[2])
                count += 1
            if count == 2:
                energies.append( energy)

                count = 0
    if idx == 1:
        energies = 0
    return np.asarray( energies )
   # except:
     #   print "error loading energies", filename_scoring
    #    return -1


def checkNone( data ):
    if data is not None:
        return data
    else:
        return 1000000


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
ext_irmsd= "-irmsd.dat"
ext_fnat = "-fnat.dat"
ext_scored = "-sorted-dr.dat"
filename_pdb = "-receptor-for-docking"

def mode_ext( num_modes):
    return "-" + num_modes +"-modes.dat"



def getCapriRanking( Irmsd, fnat):
        if fnat >= 0.5 and Irmsd <= 1.0:
            return 3
        elif fnat >= 0.3 and Irmsd <= 2.0:
            return 2
        elif fnat >= 0.1 and fnat < 0.3 and Irmsd < 4:
            return 1
        else:
            return 0

def getCapriList(list_Irmsd, list_fnat, maxIndex):
    list_capri = np.zeros(maxIndex)
    if maxIndex > len(list_Irmsd) or maxIndex  > len(list_fnat):
        print "the index is higher than the length of inpuarrays"
    else:
        for i in range(maxIndex):
            list_capri[i] = getCapriRanking(list_Irmsd[i], list_fnat[i])
    return list_capri

def capriCount(list_capri):
    count_1 = 0
    count_2 = 0
    count_3 = 0
    for rank in list_capri:
        if rank == 1:
            count_1 += 1
        elif rank == 2:
            count_2 += 1
        elif rank == 3:
            count_3 += 1
    return count_1, count_2, count_3

def evaluate( bechmarks, dict_indices):
    result = {}
    count = 0
    for name_protein, input, name_bench in bechmarks:
        count += 1
        print count,name_bench[0][0], name_protein
        result_protein = [None]*len(dict_indices)
        for name_singleBench in name_bench:
            load = True
            index_modes = name_singleBench[0].find("modes")-1
            num_modes = name_singleBench[0][index_modes]
            #print "\t",name_singleBench[0]
            if name_singleBench[0] != 'input':
                filename_rmsd = os.path.join(name_singleBench[1],name_protein+filename_pdb+ext_rmsd)
                try:
                    with open(filename_rmsd, 'r' ) as f:
                        lines = f.readlines()
                        len_list = len(lines)
                        if len_list ==0:
                            rmsd = np.zeros(1)
                            rmsd[0] = 100000
                            pos= np.zeros(1)
                        else:
                            rmsd = np.zeros( len(lines), dtype = float)
                            pos = np.arange( len(lines), dtype = int  )
                            for i, line in enumerate(lines):
                                rmsd[i] = float( line.split()[-1])
                        if rmsd is None:
                            rmsd = np.asarray(100000)
                except:
                    print "error loading rmsd file", filename_rmsd

                filename_irmsd = os.path.join(name_singleBench[1], name_protein + filename_pdb + ext_irmsd)
                irmsd = np.zeros(1)
                try:
                    with open(filename_irmsd, 'r') as f:
                        lines = f.readlines()
                        len_list = len(lines)
                        if len_list == 0:
                            irmsd = np.zeros(1)
                            irmsd[0] = 100000
                            pos = np.zeros(1)
                        else:
                            irmsd = np.zeros(len(lines), dtype=float)
                            pos = np.arange(len(lines), dtype=int)
                            for i, line in enumerate(lines):
                                irmsd[i] = float(line.split()[-1])
                        if irmsd is None:
                            irmsd = np.asarray(100000)
                except:
                    print
                    "error loading Irmsd file", filename_irmsd

                filename_fnat = os.path.join(name_singleBench[1], name_protein + filename_pdb + ext_fnat)
                fnat = np.zeros(1)
                try:
                    with open(filename_fnat, 'r') as f:
                        lines = f.readlines()
                        len_list = len(lines)
                        if len_list == 0:
                            fnat = np.zeros(1)
                            fnat[0] = 100000
                            pos = np.zeros(1)
                        else:
                            fnat = np.zeros(len(lines), dtype=float)
                            pos = np.arange(len(lines), dtype=int)
                            for i, line in enumerate(lines):
                                fnat[i] = float(line.split()[0])
                        if fnat is None:
                            fnat = np.asarray(100000)
                except:
                    print
                    "error loading Irmsd file", filename_fnat
                #plt.plot( rmsd[:20], pos[:20])
                #plt.show()
                #print "\t\t", "total rmsd", rmsd.mean(), "\n\t\t10 rmsd",rmsd[:10].mean(), "\n\t\t50 rmsd",rmsd[:50].mean(), "\n\t\t10 sorted rmsd",np.sort(rmsd)[:10].mean(),"\n\t\t50 sorted rmsd", np.sort(rmsd)[:50].mean()


                if int(num_modes) > 0:
                    filename_modes  = os.path.join(input[0],name_protein+filename_pdb+mode_ext(num_modes))
                    try:
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
                    except:
                        print "error loading mode file", filename_modes
                        amplitude = -1
                        eigenvalues = -1
                else:
                    amplitude = 0
                    eigenvalues = 0
                filename_scoring=os.path.join(name_singleBench[1],name_protein+filename_pdb+ext_scored)
                energies = getEnergyfromFile(filename_scoring)
                #print "energies",energies
            #try:

            list_capri = getCapriList(irmsd, fnat, maxIndex=50)
            one_star, two_star, three_star = capriCount(list_capri= list_capri)
            result_protein[dict_indices['mode']]    = int(num_modes)
            result_protein[dict_indices['rmsd']] = checkNone(rmsd)
            result_protein[dict_indices['irmsd']] = checkNone(irmsd)
            result_protein[dict_indices['fnat']] = checkNone(fnat)
            result_protein[dict_indices['capri']] = checkNone(list_capri)
            result_protein[dict_indices['one_star_50']] = one_star
            result_protein[dict_indices['two_star_50']] = two_star
            result_protein[dict_indices['three_star_50']] = three_star

            result_protein[dict_indices['energy']] = checkNone(energies)
            result_protein[dict_indices['pos']] = checkNone(pos)
            result_protein[dict_indices['min']] = checkNone(min(rmsd))
            result_protein[dict_indices['mean']] = checkNone(rmsd.mean())
            result_protein[dict_indices['mean_10']] = checkNone(rmsd[:10].mean())
            result_protein[dict_indices['mean_50']] = checkNone(rmsd[:50].mean())
            result_protein[dict_indices['mean_sort_10']] = checkNone(np.sort(rmsd)[:10].mean())
            result_protein[dict_indices['mean_sort_50']] = checkNone(np.sort(rmsd)[:50].mean())
            result_protein[dict_indices['amp']] = checkNone(amplitude)
            result_protein[dict_indices['ev']] = checkNone(eigenvalues)
            #except:
            #    load = False
             #   print "could not load Protein ", name_protein

           # result_protein.append( (int(num_modes), rmsd, energies,pos,  min(rmsd), rmsd.mean(), rmsd[:10].mean(), rmsd[:50].mean(), np.sort(rmsd)[:10].mean(), np.sort(rmsd)[:50].mean(),  amplitude, eigenvalues ))

        #if load:
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

    def getDataProtein(self, name_benchmark, name_protein, name_index = None):
        if name_index is None:
            return self._dict_dataBenchmark[name_benchmark][name_protein]
        else:
            return self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index[name_index]]

    #def getDataProtein(self, name_benchmark, name_protein, name_index):
    #    return self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index[name_index]]

    def setActiveBenchmarks(self, name_benchmark, name_benchmarkReference =None):
        self._activeBM = name_benchmark
        self._activeBMReference = name_benchmarkReference

    def setActiveProtein(self, name_benchmark):
        self._activeBM = name_benchmark

    def setActiveProtein(self, name_protein):
        self._activeProtein = name_protein

    def plotRmsd(self, max_indx = None, name_benchmark = None, name_protein = None, use_energy = False, show = True, clear = False):
        if name_protein is None:
            name_protein = self._activeProtein
        if name_benchmark is None:
            name_benchmark = self._activeBM
        if max_indx is None:
            max_indx = len(self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index['rmsd']])
        rmsd = self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index['rmsd']][:max_indx]
        if not use_energy:
            y = self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index['pos']][:max_indx]
        else:
            y = self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index['energy']][:max_indx]

        #plt.scatter( rmsd, y)
        if clear:
            plt.clf()
        plt.plot(rmsd, y, 'bo')
        plt.xlabel('RMSD (A)')
        plt.ylabel('Energy ')
        plt.title('{}\t{}'.format(name_benchmark, name_protein))
        if show:
            plt.show()
        return plt

    def getProteinNames(self):
        return self._proteinNames

    def getMinimum(self, name_index, name_benchmark = None, index_subarray = None, onlygetName = False):
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
                return  min( dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]]))[0]
            else:
                return  min(dict_data.items(), key=lambda x: (x[1][self._dict_index[name_index]][index_subarray]))[0]

    def getMaximum(self, name_index, name_benchmark = None, index_subarray = None, onlygetName = False):
        if name_benchmark is None:
            name_benchmark = self._activeBM
        dict_data = self._dict_dataBenchmark[name_benchmark]
        print "max ", dict_data
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

    def scoringPerformance(self, name_benchmark, name_protein):
        unsorted = self.getDataProtein(name_benchmark, name_protein, 'mean_10')
        sorted = self.getDataProtein(name_benchmark, name_protein, 'mean_sort_10')
        return sorted/ unsorted

    def getWeightedResult(self, name_benchmark, name_protein):
        try:
            data = self.getDataProtein(name_benchmark, name_protein)
            rmsdMeanBulk = np.mean(checkNone(data[self._dict_index['rmsd']][10:]))
            rmsdMeanBest = data[self._dict_index['mean_10']]
            rmsdBest = min(checkNone(data[self._dict_index['rmsd']]))

            energyMeanBulk = np.mean(checkNone(data[self._dict_index['energy']][10:]))
            energyMeanBest = np.mean(checkNone(data[self._dict_index['energy']][:10]))
            energyBest = min(checkNone(data[self._dict_index['energy']]))

            scoreRmsd = 1- (rmsdMeanBest - rmsdBest) / (rmsdMeanBulk - rmsdBest)
            scoreEnergy = 1 - (energyMeanBest - energyBest) / (energyMeanBulk - energyBest)
            #print "bulk",rmsdMeanBulk,"mean best",rmsdMeanBest,"best",rmsdBest,"energy mean bulk ",  energyMeanBulk, " energy mean best ", energyMeanBest, " energy best ", energyBest
            #print  checkNone(scoreRmsd, scoreEnergy)
            return checkNone(scoreEnergy * scoreRmsd)
        except:
            print "error getting score for benchmark ", name_benchmark, name_protein
            return 100000


    def getWeightedResultTill(self, name_benchmark, name_protein):
        try:
            data = self.getDataProtein(name_benchmark, name_protein)
            rmsdBest = min(checkNone(data[self._dict_index['rmsd']]))

            energyMeanBulk = np.mean(checkNone(data[self._dict_index['energy']]))
            energyBest = min(checkNone(data[self._dict_index['energy']]))

            siteBulk = []
            site =[]
            #print "Ess ",energyMeanBulk,  energyBest
            for E, rmsd in zip(data[self._dict_index['energy']], data[self._dict_index['rmsd']]):
                #check = 0
                s = ( E - energyMeanBulk ) / energyBest
                if rmsd < (rmsdBest + 1.5) and 0 < rmsd:
                    #print "site ", s, rmsd
                    site.append( (s,rmsd) )
                elif (rmsdBest + 1.5) < rmsd and rmsd < 1000:
                   # print "bulk ", s, rmsd
                   # if s < check:
                        #check = s
                    siteBulk.append( (s,rmsd) )
            siteMin = max(site, key=lambda t: t[0])
            siteBulkMin = max(siteBulk, key=lambda t: t[0])
           # print "check" , check
           # print siteMin, siteBulkMin
            return siteMin[0] - siteBulkMin[0]
        except:
            print "error getting score for benchmark ", name_benchmark, name_protein
            return 100000

    def stdDevRmsd(self, name_benchmark, name_protein, maxIdx):
        data = self.getDataProtein(name_benchmark, name_protein)
        rmsd = checkNone(data[self._dict_index['rmsd']])
        stdDev = 0
        mean = 0
        count = 0
        # for r in rmsd[:maxIdx]:
        #     count += 1
        #     mean += r
        # mean /= count
        # for r in rmsd[:50]:
        #     stdDev += (r - mean)* (r - mean)
        # stdDev = np.sqrt( stdDev / count)
        stdDev = np.std(rmsd[:maxIdx])
        return stdDev


    def getNumSmaller(self, name_benchmark, name_protein, borderRmsd,maxindex = None, sorted = False):
        data = self.getDataProtein(name_benchmark, name_protein)
        if maxindex is None:
            rmsd = checkNone(data[self._dict_index['rmsd']])
        else:
            rmsd = checkNone(data[self._dict_index['rmsd']][:maxindex])
        if sorted:
            rmsd = np.sort(rmsd)
        count = np.where(rmsd < borderRmsd)
        return len(count[0])

    def getNumGreater(self, name_benchmark, name_protein, borderRmsd,maxindex = None, sorted = False):
        data = self.getDataProtein(name_benchmark, name_protein)
        if maxindex is None:
            rmsd = checkNone(data[self._dict_index['rmsd']])
        else:
            rmsd = checkNone(data[self._dict_index['rmsd']][:maxindex])
        if sorted:
            rmsd = np.sort(rmsd)
        count = np.where(rmsd > borderRmsd)
        return len(count[0])
