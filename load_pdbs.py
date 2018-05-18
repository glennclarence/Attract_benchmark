import os


# This Class is used to choose categorize and identify the proteins
class ProteinEsemble:
    def __init__(self, path_Ensemble):
        self.path_ensemble = path_Ensemble
        self.list_proteins={}
        self.receptor = None
        self.ligand= None

    def add_Protein(self, filename_protein, name_protein ):
        self.list_proteins[ name_protein ] = filename_protein

    def delete_Protein(self, name_protein ):
        del self.list_proteins[ name_protein ]

    def make_Receptor(self, name_protein):
        if self.list_proteins[name_protein]:
            self.receptor =  name_protein
        for key, value in self.list_proteins.iteritems():
            if key != name_protein:
                self.ligand= key

    def make_receptor(self, filename_protein, arb):
        name_receptor = [key for key, value in self.list_proteins.iteritems() if value == filename_protein]
        if name_receptor is None:
            exit(1)
        else:
            self.receptor = name_receptor
        for key, value in self.list_proteins.iteritems():
            if value != filename_protein:
                self.ligand = key

    def get_proteinSize(self, name_protein):
        pdb = open( self.list_proteins[ name_protein])
        count_lines = len(pdb.readlines())
        return count_lines

    def get_proteinMaxSize(self):
        count_lines = 0
        maxKey = None
        for key, value in self.list_proteins.iteritems():
            pdb = open( self.list_proteins[ key ])
            length_pdb = len(pdb.readlines())
            if length_pdb > count_lines:
                count_lines = length_pdb
                maxKey  =  key
        return maxKey

    def get_pathEnsemble(self):
        return self.path_ensemble

    def get_filenameProtein(self, name_protein):
        return self.list_proteins[name_protein]

    def get_filenameAll(self):
        return [value for key, value in self.list_proteins.iteritems()]

    def get_proteinNameAll(self):
        return [key for key, value in self.list_proteins.iteritems()]

    def get_filenameReceptor(self):
        return self.list_proteins[self.receptor]

    def get_receptor(self):
        return self.receptor

    def get_ligand(self):
        return self.ligand

    def get_filenameLigand(self):
        return self.list_proteins[self.ligand]

    def get_ensemble(self):
        return self.list_proteins

def load_fromFolder( path_folder, filename_sheme):
    ensembleList = list()
    for root, dirs, files in os.walk( path_folder ):
        protein_ensemble = ProteinEsemble(root)
        if len(files)== 0:
            continue
        count = 0;
        for name in files:
            if name.endswith( filename_sheme ):
                count +=1
                filename_protein = os.path.join( root, name )
                protein_ensemble.add_Protein( filename_protein, name.split('.')[0] )
        if count > 0:
            ensembleList.append( protein_ensemble )
    return ensembleList

def set_receptorBySize( protein_ensemble ):
    receptorName = protein_ensemble.get_proteinMaxSize()
    protein_ensemble.make_Receptor(receptorName)


