import os


# This Class is used to choose categorize and identify the proteins
class ProteinEsemble:
    def __init__(self, path_Ensemble):
        self.path_ensemble = path_Ensemble
        self.list_PDB ={}
        self.receptor = None
        self.ligand= None

    def add_Protein(self,  name_protein, filename_PDB ):
        self.list_PDB[ name_protein ] = filename_PDB

    def delete_Protein(self, name_protein ):
        del self.list_PDB[name_protein]

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

    def get_filenamePDB(self, name_protein):
        return self.list_PDB[name_protein]

    def get_ensemblePDB(self):
        return self.list_PDB


def load_fromFolder( path_folder, filename_sheme_pdb  ):
    ensembleList = list()
    for root, dirs, files in os.walk( path_folder ):
        protein_ensemble = ProteinEsemble(root)
        if len(files)== 0:
            continue
        count = 0;
        for name in files:
            if name.endswith( filename_sheme_pdb ) and filename_sheme_pdb is not None:
                count += 1
                filename_pdb = os.path.join( root, name )

                name_protein = os.path.split(root)[-1] + "-" + name.split('.')[0]
                protein_ensemble.add_Protein(name_protein, filename_pdb)

        if count > 0:
            ensembleList.append( protein_ensemble )
    return ensembleList

def get_receptorBySize( protein_ensemble ):
    receptorName = get_proteinMaxSize(protein_ensemble)
    return receptorName

def get_ligandBySize(protein_ensemble):
    ligandName = get_proteinMinSize(protein_ensemble)
    return ligandName

def get_proteinSize(self, name_protein):
    pdb = open( self.list_proteins[ name_protein])
    count_lines = len(pdb.readlines())
    return count_lines

def get_proteinMaxSize( protein_list):
    count_lines = 0
    maxKey = None
    for key, value in protein_list.iteritems():
        pdb = open( protein_list[ key ])
        length_pdb = len(pdb.readlines())
        if length_pdb > count_lines:
            count_lines = length_pdb
            maxKey  =  key
    return maxKey

def get_proteinMinSize( protein_list):
    count_lines = 10000000
    maxKey = None
    for key, value in protein_list.iteritems():
        pdb = open( protein_list[ key ])
        length_pdb = len(pdb.readlines())
        if length_pdb < count_lines:
            count_lines = length_pdb
            maxKey  =  key
    return maxKey