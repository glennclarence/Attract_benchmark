



class ResultClass:
    def __int__(self, dict_index):

        self._dict_index = dict_index
        self._dict_dataBenchmark = {}
        self._activeBM = None
        self._activeBMReference = None
        self._activeProtein = None
        self._proteinNames = set()

    def add_dataBenchmark(self, name_benchmark, data_benchmark):
        self._dict_dataBenchmark[name_benchmark] = data_benchmark
        for name_protein, data in data_benchmark:
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

    def plotRmsd(self,name_benchmark = None, name_protein = None, use_energy = False):
        if name_protein is None:
            name_protein = self._activeProtein
        if name_benchmark is None:
            name_benchmark = self._activeBM
        rmsd = self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index['idx_rmsd']]
        if not use_energy:
            y = self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index['idx_pos']]
        else:
            y = self._dict_dataBenchmark[name_benchmark][name_protein][self._dict_index['idx_energy']]
        plt.plot( rmsd, y)
        plt.show()
    def getProteinNames(self):
        return self._proteinNames

    def getMininum(self, name_index, name_benchmark = None):

        if name_benchmark is None:
            name_benchmark = self._activeBM
        dict_data = self._dict_dataBenchmark[name_benchmark]
        return min(dict_data[name_benchmark], key=dict_data[name_benchmark].get)

    def getMaximum(self, name_index, name_benchmark = None):

        if name_benchmark is None:
            name_benchmark = self._activeBM
        dict_data = self._dict_dataBenchmark[name_benchmark]
        return max(dict_data[name_benchmark], key=dict_data[name_benchmark].get)