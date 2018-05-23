import  load_pdbs as lpdb
import ProteinConfiguration as protconf
import os
import run_computation as compute
import run_analysis as analyse
from run_analysis import join_modefiles
import benchmark_timer as benchtime
import time

class ProteinPair:
    def __init__(self, receptor, ligand, filename_dof, input_folder, output_folder, filename_modesJoined = None ):
        self.receptor = receptor
        self.ligand = ligand
        self.filename_dof =filename_dof
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.filename_modesJoined = filename_modesJoined
    def get_receptor(self):
        return self.receptor
    def get_ligand(self):
        return self.ligand
    def get_filenameDOF(self):
        return self.filenameDOF
    def get_filenameModes(self):
        return self.filename_modesJoined


def configure_protein( path_inputFolder, path_outputFolder, name_protein, filename_pdb_protein, chain = None, num_modes = 0):
    protein = protconf.ProteinConfiguration(filename_pdb_protein= filename_pdb_protein, name_protein=name_protein)
    protein.set_path_inputFolder( path_inputFolder = path_inputFolder )
    protein.set_path_outputFolder( path_outputFolder = path_outputFolder )
    protein.set_chain( chain )
    protein.set_num_modes( num_modes = num_modes )
    return protein

# receptor = protconf.ProteinConfiguration( filename_pdb_protein=filename_receptor, name_protein=name_receptor)
#         receptor.set_path_inputFolder( path_inputFolder=os.path.join(ensemble.get_pathEnsemble(), 'input'))
#         receptor.set_path_outputFolder( path_outputFolder=os.path.join(ensemble.get_pathEnsemble(), name_benchmark))
#         receptor.set_filename( filename_reduce= filename_receptor)
#         receptor.set_chain(chain_receptor)
#         receptor.set_num_modes(num_modes)
# ligand = protconf.ProteinConfiguration( filename_pdb_protein=filename_ligand, name_protein=name_ligand)
#         ligand.set_path_inputFolder( path_inputFolder=os.path.join(ensemble.get_pathEnsemble(), 'input'))
#         ligand.set_path_outputFolder( path_outputFolder=os.path.join(ensemble.get_pathEnsemble(), name_benchmark))
#         ligand.set_filename( filename_reduce=filename_ligand)
#         ligand.set_chain(chain_ligand)
#         ligand.set_num_modes(num_modes)



# running a benchmark is divided in four steps:

#                   Identifying -------------> configuring -------------> computing -------------> analysis
#
# 1.step: use load_fromFolder() from load_pdbs to indentify all the pdbs that match a scheme. This function returns a list of type ProteinEsemble 's containing the folder and the filename of the receptor and the ligand
# 2.step: by using the filenames of te proteinensembles you can create a protein configuration for each protein.
# By doing that for the ligand , the receptor and by creating a file with the starting dof with ProteinConfiguration.create_StartingPositions() you can create a pair.
#3.step: create an object of type run_computation.Worker for the docking and the scoring step and add all the pairs to the worker via worker.add_ensembleToQueue() and start computation via worker.comput_serial().
# There is multithreaded version available also using worker.start_workers() and worker.run()
# 4. step: run analysis via run_analysis.run_analysis()




def run_benchmark( path_folder, filename_scheme, name_benchmark, create_grid = False, create_modes = False, create_dofs = False, create_reduce =False, num_modes= 0, use_orig = False, num_threads = 7 ):
    pairs = list()

    benchmark = benchtime.measure_benchmark()
    benchmark.timer_add("Create_reducepdb")
    benchmark.timer_add("Create_modes")
    benchmark.timer_add("Create_grid")
    benchmark.timer_add("Create_dofs")
    benchmark.timer_add("Minimization")
    benchmark.timer_add("Scoring")
    benchmark.timer_add("Analysis")

    pdb = False
    pdb_reduced = True
    pdb_allatom = False
    analyse = False
    do_scoring = False
    index_chain = 4

    parameter_dir = os.environ['ATTRACTDIR'] + "/../attract.par"
    filename_extension_dock = "_dock.result"
    filename_extension_scoring = "_scoring.result"

    #dock = compute.Worker(path_attract="/home/glenn/Documents/Masterarbeit/git/gpuATTRACT_2.0", do_minimization=True, num_threads=num_threads , args = benchmark, use_OrigAttract= use_orig)
    #score = compute.Worker(path_attract="/home/glenn/Documents/Masterarbeit/git/gpuATTRACT_2.0", do_scoring=True, num_threads= num_threads , args = benchmark, use_OrigAttract= use_orig)

    dock = compute.Worker(path_attract="/home/glenn/Documents/attract/bin", name_attractBinary = "attract", do_minimization=True,
                          num_threads=num_threads, args=None, use_OrigAttract=use_orig)
    score = compute.Worker(path_attract="/home/glenn/Documents/attract/bin", name_attractBinary = "attract", do_scoring=True,
                           num_threads=num_threads, args=None, use_OrigAttract=use_orig)


    print '**************************************************************'
    print "Load Protein Pdbs"
    #load all proteins inside one folder into
    protein_ensembles =         lpdb.load_fromFolder( path_folder= path_folder, filename_sheme_pdb=filename_scheme)
    protein_ensembles_allatom = lpdb.load_fromFolder( path_folder=path_folder, filename_sheme_pdb="-aa.pdb")

    print 'detected',len(protein_ensembles) , "protein ensembles"
    print '**************************************************************'
    print "Create Protein Configurations"
    for count, ensemble in enumerate(protein_ensembles):
        
        proteins = ensemble.get_ensemblePDB()
        path_input = os.path.join(ensemble.get_pathEnsemble(), 'input')
        path_output = os.path.join(ensemble.get_pathEnsemble(), name_benchmark)
        filename_dof = os.path.join(os.path.join(ensemble.get_pathEnsemble(), 'input'), "dofs.dat")

        name_receptor = lpdb.get_receptorBySize(   proteins )
        filename_receptor = proteins[ name_receptor]
        chain_receptor = protconf.get_chainfromName(name_receptor, index_chain)
        receptor = configure_protein(path_inputFolder = path_input,     path_outputFolder = path_output, name_protein = name_receptor,
                                     filename_pdb_protein = None,       chain = chain_receptor,          num_modes = num_modes)
        receptor.set_filename(       filename_reduce= filename_receptor )


        name_ligand = lpdb.get_ligandBySize(proteins)
        filename_ligand = proteins[name_ligand]
        chain_ligand = protconf.get_chainfromName(name_ligand, index_chain)
        ligand = configure_protein( path_inputFolder=path_input,    path_outputFolder=path_output,  name_protein=name_ligand,
                                    filename_pdb_protein=None,      chain=chain_ligand,             num_modes=num_modes)
        ligand.set_filename(            filename_reduce=filename_ligand)


        print "{}/{} ".format( count, len(protein_ensembles)),"\t--Configure protein: ", name_receptor, " and ", name_ligand
        if create_reduce:
            print "\t\t--Create reduced pdbs"
            benchmark.timer_start("Create_reducepdb")
            ligand.reduce( overwrite=False)
            receptor.reduce( overwrite=False )
            benchmark.timer_appendStop("Create_reducepdb")
        if create_grid:
            print "\t\t--Create grids"
            ligand.set_partner( receptor.get_filenamePdbReduced() )
            receptor.set_partner(ligand.get_filenamePdbReduced())
            benchmark.timer_start("Create_grid")
            receptor.create_grid(overwrite=False)
            ligand.create_grid(overwrite=False)
            benchmark.timer_appendStop("Create_grid")
        if create_modes:
            print "\t\t--Create modefiles"
            benchmark.timer_start("Create_modes")
            receptor.create_modes(overwrite=False)
            ligand.create_modes(overwrite=False)
            benchmark.timer_appendStop("Create_modes")
            filename_modesJoined = None
            if use_orig:
                filename_modesJoined = os.path.join( receptor.get_pathInput(), "allModes.dat" )
                join_modefiles( receptor.get_filenameModes(), ligand.get_filenameModes(), filename_modesJoined )



        if create_dofs:
            print "\t\t--Create dofs"
            benchmark.timer_start("Create_dofs")
            protconf.create_StartingPositions("/home/glenn/Documents/attract/rotation.dat",  receptor.get_filenamePdbReduced(), ligand.get_filenamePdbReduced(),filename_dof)
            benchmark.timer_appendStop("Create_dofs")

        pair = ProteinPair(receptor, ligand, filename_dof, receptor.get_pathInput(), receptor.get_pathOutput() , filename_modesJoined = filename_modesJoined )
        pairs.append( pair )



    parameter_dir = os.environ['ATTRACTDIR'] + "/../attract.par"

    print '**************************************************************'
    print "Run Minimization"
    for count, pair in enumerate(pairs):
        receptor = pair.get_receptor()
        ligand = pair.get_ligand()
        print "{}/{} ".format( count, len(pairs)),"\t--minimize protein: ", receptor.get_name(), " and ", ligand.get_name()
        filename_output = os.path.join(receptor.get_pathOutput(), receptor.get_name() + filename_extension_dock)

        dock.add_ensembleToQueue(  filename_dofs=pair.filename_dof, filename_output=filename_output,    filename_parameter=parameter_dir,
                                  filename_pdbReceptor=receptor.get_filenamePdbReduced(),               filename_alphabetReceptor=receptor.get_filenameAlphabet(),
                                  filename_gridReceptor=receptor.get_filenameGrid(),                    filename_modesReceptor=receptor.get_filenameModes(),
                                  num_modesReceptor=num_modes,      filename_pdbLigand=ligand.get_filenamePdbReduced(), filename_alphabetLigand=ligand.get_filenameAlphabet(),
                                  filename_gridLigand=ligand.get_filenameGrid(), filename_modesLigand=ligand.get_filenameModes(),
                                  num_modesLigand=num_modes, filename_modesJoined= pair.get_filenameModes() )

    benchmark.timer_start("Minimization")
    dock.start_threads()
    benchmark.timer_appendStop("Minimization")



    if do_scoring:
        print '**************************************************************'
        print "Run Scoring"
        for count, pair in enumerate(pairs):
            receptor = pair.get_receptor()
            ligand = pair.get_ligand()
            print "{}/{} ".format( count, len(pairs)),"\t--score protein: ", receptor.get_name(), " and ", ligand.get_name()
            filename_output = os.path.join(receptor.get_pathOutput(), receptor.get_name() + filename_extension_scoring)

            filename_docking = os.path.join(pair.receptor.get_pathOutput(),
                                            pair.receptor.get_name() + filename_extension_dock)

            score.add_ensembleToQueue(  filename_dofs=filename_docking,filename_output=filename_output,filename_parameter=parameter_dir,
                                      filename_pdbReceptor=receptor.get_filenamePdbReduced(),filename_alphabetReceptor=receptor.get_filenameAlphabet(),
                                      filename_gridReceptor=receptor.get_filenameGrid(),filename_modesReceptor=receptor.get_filenameModes(),
                                      num_modesReceptor=num_modes,filename_pdbLigand=ligand.get_filenamePdbReduced(),filename_alphabetLigand=ligand.get_filenameAlphabet(),
                                      filename_gridLigand=ligand.get_filenameGrid(),filename_modesLigand=ligand.get_filenameModes(),
                                      num_modesLigand=num_modes, filename_modesJoined= pair.get_filenameModes() )

    benchmark.timer_start("Scoring")
    benchmark.timer_appendStop("Scoring")

    score.stop_threads_if_done()
    dock.stop_threads_if_done()
    print '**************************************************************'
    print "Run Analysis"



    if analyse is True:
        for count, pair in enumerate( pairs ):
            allatom_proteins = protein_ensembles_allatom[count].get_ensemblePDB()
            filename_scoring = os.path.join(pair.receptor.get_pathOutput(), pair.receptor.get_name() + filename_extension_scoring)
            filename_docking = os.path.join(pair.receptor.get_pathOutput(), pair.receptor.get_name() + filename_extension_dock)

            print "{}/{} ".format( count, len(pairs)),"\t--analyse protein: ", pair.receptor.get_name(), " and ", pair.ligand.get_name()
            #allatom_proteins["receptor-aa"],allatom_proteins["ligand-aa"]
            analyse.run_analysis(pair.output_folder,pair.receptor.name_protein, filename_docking, filename_scoring,  pair.receptor.get_filenamePdbReduced(),pair.ligand.get_filenamePdbReduced(),pair.ligand.get_filenamePdbReduced(),filename_modesReceptor =pair.receptor.get_filenameModes(),
                                 filename_modesLigand=pair.ligand.get_filenameModes(),
                     num_modesReceptor=num_modes, num_modesLigand=num_modes, filename_modesJoined = filename_modesJoined,
                     path_attract=os.environ['ATTRACTDIR'], path_attractTools=os.environ['ATTRACTTOOLS'])

    benchmark.save_benchmark( os.path.join(path_folder, name_benchmark + '-time.dat'))

run_benchmark( "/home/glenn/Documents/Attract_benchmark2", "-unboundr.pdb",name_benchmark = "benchmark_5_modes_orig", create_grid = True, create_modes = True, create_dofs = True, create_reduce = False, num_modes = 5, use_orig= True)

