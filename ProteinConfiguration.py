import os, sys
import subprocess
sys.path.append('/home/glenn/Documents/Masterthesis/attract_unchanged/allatom')
#import aareduce_01.py

#sys.path.append('/home/glenn/Documents/Masterthesis/attract_unchanged/tools/custom/')

from reduce_nonparse import reduce
from modes_nonparse import modes
from aareduce_nonparse import aareduce

#This filecontains two classes and one free function:

#1. completePath is a simple class which offers safe handling of files and paths.
# The class holds a path and filnames can either be added to the path of is can be checked if a file exists.

#2. ProteinConfiguration : This class deals with the complete configuration of a single Protein prior to docking.
#   1. setting the paths that are needed such as input and outputpath
#   2. setting the filenames based on the name of the protein/configuration for the reduction modefiles, grids and so on
#   3. creating all the files that are necessary to dock the protein such as the modefile, gridfile, alphabet file and reduced Pdb

#3. create_StartingPositions: creates the starting degrees of freedom for two protein pdbs


class completePath:
    def __init__(self, path):
        self.path = path
    def filename(self, filename):
        return os.path.join( self.path, filename )
    def isfile(self, filename):
        return os.path.isfile( os.path.join( self.path, filename))



#to init
class ProteinConfiguration:
    """Sets the configuration of a Protein according to a name(also the beginning of the created files) and the filename of the pdbfile.
    The configuration includes all the files that are necessary for docking such as reduced pdb, modesfile, gridfile and alphabet file.
    To complete the configuration the following methods can be called:
    1. set_filename: allows to set names of the configuration files individually
    2. set_path_inputFolder: sets and creates a path where all the configuration files are stored in
    3. set_path_outFolder: sets and creates a path where all the results of the docking are stored in
    4. set_num_modes: sets the number of modes( default:0 -> no modefile is created)
    5. set_path_attract: path to binaries such as systsearch, translate and make-grid-omp (default: attract/bin) and the parameter file attract.par
    6. set_chain: sets the chain a protein belongs to ( default = "A")
    7. set_partner: sets the partner for the protein in form a reduced pdb file, which is need later on for the creation of the grid.

    After setting all paths and filenames, files can be created via reduce, create_modes and create_grid
    """
    def __init__(self, name_protein, filename_pdb_protein, path_attract = os.environ['ATTRACTDIR'] , path_attractTools = os.environ['ATTRACTTOOLS']):
        self.name_protein = name_protein
        self.use_modes = False

        self.path_attract = path_attract
        self.path_attractTools = path_attractTools
        self.ext_mapping = ".mapping"
        self.ext_reduce = "_reduce.pdb"
        self.ext_aa = "_aa.pdb"
        self.ext_modes = "_modes.dat"
        self.ext_alphabet = "_grid.alphabet"
        self.ext_grid = "grid.grid"
        self.chain = "A"

        self.filename_pdb_protein = filename_pdb_protein
        self.filename_reduce    = name_protein + self.ext_reduce
        self.filename_alphabet  = name_protein + self.ext_alphabet
        self.filename_grid      = name_protein + self.ext_grid
        self.filename_modes     = name_protein + self.ext_modes
        self.filename_allAtom   = name_protein + self.ext_aa





    def set_filename( self, filename_reduce = None, filename_modes = None, filename_allAtom = None,
                      filename_alphabet=None, filename_grid = None ):
        """allows to set the filenames individually. By default the filenames are set to name_protein + self.extension.
        The filename can either include the whole path or just the filename ( in case it is located in the inputFolder"""
        if filename_reduce is not None:
            self.filename_reduce = filename_reduce
        if filename_alphabet is not None:
            self.filename_alphabet = filename_alphabet
        if filename_grid is not None:
            self.filename_grid = filename_grid
        if filename_modes is not None:
            self.filename_modes = filename_modes
        if filename_allAtom is not None:
            self.filename_allAtom = filename_allAtom
    #not used right now
    def create_file_from_list( self, filename, list ):
        """ writes the key and its value to a file"""
        content = ""
        for item in list:
            content += "{}\n".format( item )

        file = open( filename, "w")
        file.write(content)
        file.close()

    def get_name(self):
        """returns the name of the configuration."""
        return self.name_protein

    def set_path_inputFolder( self, path_inputFolder ):
        """sets and creates the folder where all the configuration files are created in."""
        self.path_inputFolder  = path_inputFolder
        bash_command = "mkdir -p {} ".format( self.path_inputFolder )
        os.system( bash_command )

    def set_path_outputFolder( self, path_outputFolder ):
        """sets and creates the folder where all the result files are created in."""
        self.path_outputFolder  = path_outputFolder
        bash_command = "mkdir -p {} ".format(self.path_outputFolder)
        os.system(bash_command)

    def set_num_modes( self, num_modes):
        """sets the number of modes and wheter modefiles can be created of not"""
        self.num_modes = num_modes
        self.use_modes = self.num_modes > 0

    def get_numModes(self):
        return self.num_modes

    def set_chain(self, chain):
        """sets the chain the protein belongs to. The default is 'A' which is set by __init__ ."""
        self.chain = chain
    def get_chain(self):
        return self.chain

    def set_path_attract( self,  path_attract = os.environ['ATTRACTDIR']):
        """Sets the folder to binaries that are needed for the creation of the grid  and the parameters used in the reduce method."""
        self.path_attract = path_attract

    def set_attractTools(self , path_attractTools = os.environ['ATTRACTTOOLS']):
        self.path_attractTools = path_attractTools

    def set_device( self, device= 'CPU'):
        """Either CPU can to used. Can be used to pass the computation device to the actual computation."""
        self.device = device

    def get_pathOutput(self):
        return self.path_outputFolder

    def get_pathInput(self):
        return self.path_inputFolder

    def get_filenamePdbReduced(self):
        return os.path.join( self.path_inputFolder, self.filename_reduce)

    def get_filenameGrid(self):
        return os.path.join( self.path_inputFolder, self.filename_grid)

    def get_filenameAlphabet(self):
        return os.path.join( self.path_inputFolder, self.filename_alphabet)

    def get_filenameModes(self):
        return os.path.join( self.path_inputFolder, self.filename_modes )


    def reduce( self, overwrite = False):
        """Creates the allatom.pdb and the reduced pdb. Additionally you can choose to overwrite the files, in case they already exist. """
        path = completePath( self.path_inputFolder )
        output_path         = self.path_inputFolder+"/"
        output_name_mapping =self.name_protein+ self.ext_mapping


        if path.isfile(  self.filename_allAtom ) is False or (path.isfile(  self.filename_allAtom ) and overwrite is True):
            aareduce(args_pdb=self.filename_pdb_protein,
                output_path         = output_path,
                output_name_pdb     =  self.filename_allAtom,
                output_name_mapping = output_name_mapping,
                args_chain=self.chain,  args_pdb2pqr=True, args_dumppatch=True);
        if path.isfile(  self.filename_reduce ) is False or( path.isfile(  self.filename_reduce ) and overwrite is True):
            reduce(pdb=self.filename_pdb_protein,
                   path_output= output_path,
                   name_output =  self.filename_reduce,  chain=self.chain);

        return   path.filename(self.filename_allAtom),   path.filename(self.filename_reduce)

    def create_modes(self, overwrite = False):
        """Creates the modefile in the inputFolder according to the number of modes set before."""
        path = completePath(self.path_inputFolder)
        path_output = self.path_inputFolder + "/"
        if path.isfile(self.filename_modes) is False or ( path.isfile(self.filename_modes) and overwrite is True):
            modes(        pdb = self.get_filenamePdbReduced(),
                  path_output = path_output,
                  name_output = self.filename_modes, aapdb=None, nrmodes=5, sugar=False, scale=1.0, aacontact=False)
        return path.filename(self.filename_modes)


    def set_partner( self, filename_pdb_partner):
        """Sets a partner protein which is later on docked to. This is needed for the creation of the grid."""
        self.filename_pdb_partner = filename_pdb_partner

    def create_grid( self, overwrite = False ):
        path = completePath(self.path_inputFolder)

        bash_command = "awk '{print substr($0,58,2)}'"
        bash_command += " {} ".format( self.filename_pdb_partner )
        bash_command += " | sort -nu > {}".format(os.path.join(self.path_inputFolder, self.filename_alphabet))
        if path.isfile(self.filename_alphabet ) is False or (path.isfile(self.filename_alphabet ) and overwrite is True):
            os.system(bash_command )

        bash_command = "{}/make-grid-omp {} {}/../attract.par".format( self.path_attract, os.path.join(self.path_inputFolder, self.filename_reduce), self.path_attract);
        bash_command += " 10.0 12.0 {} --alphabet {} 2>/dev/null".format( os.path.join(self.path_inputFolder, self.filename_grid), os.path.join(self.path_inputFolder, self.filename_alphabet) );
        if path.isfile(self.filename_grid) is False or (path.isfile(self.filename_grid) and overwrite is True):
            os.system(bash_command)
        return  path.filename(self.filename_grid),  path.filename(self.filename_alphabet)


def create_StartingPositions( file_input_rotation, file_input_pdbReducedReceptor, file_input_pdbReducedLigand, file_output_DOFs, overwrite = False):
    """Creates the starting positions according to the coarse grain pdbs."""

    path_attract = os.environ[ 'ATTRACTDIR']
    path_output = os.path.split(file_output_DOFs)[0]


    if os.path.isfile(file_output_DOFs ) is False or (os.path.isfile(file_output_DOFs ) and overwrite is True):
        bash_command = "cp {} {}/{}".format(file_input_rotation, path_output, "rotation.dat")
        os.system( bash_command )
        bash_command = "{}/translate {} {} > {}/translate.dat".format( path_attract, file_input_pdbReducedReceptor, file_input_pdbReducedLigand, path_output)
        os.system(bash_command)
        bash_command = "{}/tools/systsearch {}/rotation.dat {}/translate.dat> {}".format( os.path.dirname(os.path.abspath(__file__)),path_output,path_output, file_output_DOFs )
        os.system( bash_command )


def get_chainfromName(name, index):
    """Return the chain of as a letter from a filename at index 'index'"""
    return name[index]


# not used right now
import contextlib
@contextlib.contextmanager
def redirect_argv(num):
    sys._argv = sys.argv[:]
    sys.argv=[str(num)]
    yield
    sys.argv = sys._argv

with redirect_argv(1):
    print(sys.argv)
