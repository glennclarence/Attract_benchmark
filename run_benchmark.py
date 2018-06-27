import  load_pdbs as lpdb
import ProteinConfiguration as protconf
import os
import run_computation as compute
import run_analysis as analyse
from run_analysis import join_modefiles
import benchmark_timer as benchtime
import time
from Queue import Queue

class ProteinPair:
    def __init__(self,id, receptor, ligand, filename_dof, input_folder, output_folder, filename_modesJoined = None, filename_modesJoined_reference=None, filename_pdb_reference_ligand = None ):
        self.receptor = receptor
        self.ligand = ligand
        self.filename_dof =filename_dof
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.filename_pdb_reference_ligand = filename_pdb_reference_ligand
        self.filename_modesJoined = filename_modesJoined
        self.filename_modesJoined_reference = filename_modesJoined_reference
        self.id = id
    def get_receptor(self):
        return self.receptor
    def get_ligand(self):
        return self.ligand
    def get_filenameDOF(self):
        return self.filenameDOF
    def get_filenameModes(self):
        return self.filename_modesJoined


def configure_protein( path_inputFolder, path_outputFolder, name_protein, filename_pdb_protein,filename_reference = None, chain = None, num_modes = 0):
    protein = protconf.ProteinConfiguration(filename_pdb_protein= filename_pdb_protein, name_protein=name_protein, filename_pdb_reference=filename_reference)
    protein.set_path_inputFolder( path_inputFolder = path_inputFolder )

    protein.set_path_outputFolder( path_outputFolder = path_outputFolder )
    protein.set_chain( chain )
    protein.set_num_modes( num_modes = num_modes )
    return protein

# running a benchmark is divided in four steps:

#                   Identifying -------------> configuring -------------> computing -------------> analysis
#
# 1.step: use load_fromFolder() from load_pdbs to indentify all the pdbs that match a scheme. This function returns a list of type ProteinEsemble 's containing the folder and the filename of the receptor and the ligand
# 2.step: by using the filenames of te proteinensembles you can create a protein configuration for each protein.
# By doing that for the ligand , the receptor and by creating a file with the starting dof with ProteinConfiguration.create_StartingPositions() you can create a pair.
#3.step: create an object of type run_computation.Worker for the docking and the scoring step and add all the pairs to the worker via worker.add_ensembleToQueue() and start computation via worker.comput_serial().
# There is multithreaded version available also using worker.start_workers() and worker.run()
# 4. step: run analysis via run_analysis.run_analysis()




def run_benchmark( path_folder, filename_scheme, name_benchmark, create_grid = False, create_modes = False, create_dofs = False, create_reduce =False, num_modes= 0, use_orig = False, num_threads = 7, do_analyse = True, do_scoring = True, do_minimization = True ):
    pairs = {}
    finisheditems_scoring = Queue()
    finisheditems_docking = Queue()

    benchmark = benchtime.measure_benchmark()
    benchmark.timer_add("Create_reducepdb")
    benchmark.timer_add("Create_modes")
    benchmark.timer_add("Create_grid")
    benchmark.timer_add("Create_dofs")
    benchmark.timer_add("Analysis")

    benchmark_minimization = benchtime.measure_benchmark()
    benchmark_scoring = benchtime.measure_benchmark()


    pdb = False
    pdb_reduced = True
    pdb_allatom = False

    do_configuration = True
    index_chain = 5

    parameter_dir = os.environ['ATTRACTDIR'] + "/../attract.par"
    filename_extension_dock = "_dock.result"
    filename_extension_scoring = "_scoring.result"

    path_attract_app = "/home/glenn/Documents/Masterarbeit/git/gpuATTRACT_2.0"
   # path_attract_app = "/home/glenn/Documents/Masterarbeit/git/gpuattract_backup/gpuATTRACT_2.0"
    name_attractBinary = "AttractServer"
    if use_orig:
        path_attract_app = "/home/glenn/Documents/attract/bin"

        name_attractBinary = "attract"
    #if do_minimization:
    dock = compute.Worker(path_attract=path_attract_app, name_attractBinary = name_attractBinary, do_minimization=True,
                          num_threads=num_threads, args=(benchmark_minimization,finisheditems_docking,), use_OrigAttract=use_orig)

    path_attract_app = "/home/glenn/Documents/attract/bin"
    name_attractBinary = "attract"
    dock.start_threads()
    #if do_scoring:
    score = compute.Worker(path_attract=path_attract_app, name_attractBinary= name_attractBinary,
                               do_scoring=True,
                               num_threads=num_threads, args=(benchmark_scoring, finisheditems_scoring,),
                               use_OrigAttract=True)

    score.start_threads()


    print '**************************************************************'
    print "Load Protein Pdbs"
    #load all proteins inside one folder into
    protein_ensembles =         lpdb.load_fromFolder( path_folder= path_folder, filename_sheme_pdb=filename_scheme)
    protein_ensembles_reference = lpdb.load_fromFolder( path_folder=path_folder, filename_sheme_pdb="-aa.pdb")

    print 'detected',len(protein_ensembles) , "protein ensembles"
    print '**************************************************************'
    print "Create Protein Configurations"
    for count, ensemble in enumerate(protein_ensembles):
        #set paths and get proteins from the ensemble list
        proteins = ensemble.get_ensemblePDB()
        path_input =    os.path.join( ensemble.get_pathEnsemble(), 'input')
        path_output =   os.path.join( ensemble.get_pathEnsemble(), name_benchmark)
        filename_dof =  os.path.join( os.path.join(ensemble.get_pathEnsemble(), 'input'), "dofs.dat")
        name_pair = os.path.split( ensemble.get_pathEnsemble() )[-1]
        filename_modesJoined = None



        #get the receptor and the ligand via their size and create Proteinconfigurations from it for receptor and ligand
       # name_receptor = lpdb.get_receptorBySize(   proteins )
        name_receptor = lpdb.get_proteinByindex(proteins,5, 'r')
        #name_receptor = lpdb.get_proteinByindex(proteins, 14, 'r')
        #filename_receptor_reference = os.path.join( os.path.dirname( proteins[ name_receptor] ), '1AVX_r_b.pdb' )
        filename_receptor_reference = os.path.join(os.path.dirname(proteins[name_receptor]), 'receptor-refe.pdb')
        #filename_receptor_reference = os.path.join(os.path.dirname(proteins[name_receptor]), '1AVX_r_b.pdb')
        filename_receptor = proteins[ name_receptor]

        chain_receptor = protconf.get_chainfromName(name_receptor, index_chain)
        chain_receptor = "A"
        receptor = configure_protein(path_inputFolder = path_input,     path_outputFolder = path_output, name_protein = name_receptor,
                                     filename_pdb_protein = filename_receptor,       chain = chain_receptor,          num_modes = num_modes, filename_reference=filename_receptor)
        #receptor.set_filename(       filename_reduce= filename_receptor )


        #name_ligand = lpdb.get_ligandBySize(proteins)
        name_ligand = lpdb.get_proteinByindex(proteins, 5, 'l')
        #name_ligand = lpdb.get_proteinByindex(proteins, 14, 'l')
        filename_ligand = proteins[name_ligand]
        filename_ligand_reference = os.path.join( os.path.dirname( proteins[ name_ligand] ), 'ligand-refe.pdb' )
        #filename_ligand_reference = os.path.join(os.path.dirname(proteins[name_ligand]), '1AVX_l_b.pdb')
        #filename_ligand_reference = filename_ligand
        chain_ligand = protconf.get_chainfromName(name_ligand, index_chain)
        chain_ligand = "B"
        ligand = configure_protein( path_inputFolder=path_input,    path_outputFolder=path_output,  name_protein=name_ligand,
                                    filename_pdb_protein=filename_ligand,      chain=chain_ligand,             num_modes=num_modes, filename_reference=filename_ligand)
        #ligand.set_filename(            filename_reduce=filename_ligand)

        if do_configuration:
            print "{}/{} ".format( count + 1, len(protein_ensembles)),"\t--Configure protein: ", name_receptor, " and ", name_ligand
            if create_reduce:
                print "\t\t--Create reduced pdbs"
                benchmark.timer_start("Create_reducepdb")
                ligand.reduce( overwrite=False, allatom=True)
                receptor.reduce( overwrite=False , allatom=True)
                benchmark.timer_appendStop("Create_reducepdb")
            if create_grid:
                print "\t\t--Create grids"
                ligand.set_partner( receptor.get_filenamePdbReduced() )
                receptor.set_partner(ligand.get_filenamePdbReduced())
                benchmark.timer_start("Create_grid")
                receptor.create_grid(overwrite=False)
                ligand.create_grid(overwrite=False)
                benchmark.timer_appendStop("Create_grid")
            filename_modesJoinedReference = None
            if create_modes and num_modes > 0:
                print "\t\t--Create modefiles"
                benchmark.timer_start("Create_modes")
                receptor.create_modes(overwrite=False)
                ligand.create_modes(overwrite=False)
                benchmark.timer_appendStop("Create_modes")
                filename_modesJoined = None


                filename_modesJoinedReference = os.path.join(receptor.get_pathInput(), "allModes-reference-"+str(num_modes)+".dat")
                #if use_orig:
                filename_modesJoined = os.path.join( receptor.get_pathInput(), "-allModes-"+str(num_modes)+".dat" )
                join_modefiles( receptor.get_filenameModes(), ligand.get_filenameModes(), filename_modesJoined )

                if not os.path.isfile(filename_modesJoinedReference):
                    protconf.create_modes(os.path.join(receptor.get_pathInput(),receptor.filename_allAtom), receptor.get_pathInput(), "/receptor-Modes-reference.dat", num_modes)
                    protconf.create_modes(os.path.join(ligand.get_pathInput(),ligand.filename_allAtom), ligand.get_pathInput(), "/ligand-Modes-reference.dat", num_modes)
                    #protconf.create_modes(os.path.join(receptor.get_pathInput(), filename_receptor_reference),
                    #                      receptor.get_pathInput(), "/receptor-Modes-reference.dat", num_modes)
                    #protconf.create_modes(os.path.join(ligand.get_pathInput(), filename_ligand_reference),
                    #                      ligand.get_pathInput(), "/ligand-Modes-reference.dat", num_modes)
                    join_modefiles( os.path.join(receptor.get_pathInput(),"receptor-Modes-reference.dat"), os.path.join(ligand.get_pathInput(),"ligand-Modes-reference.dat"), filename_modesJoinedReference)
            if create_dofs:
                print "\t\t--Create dofs"
                benchmark.timer_start("Create_dofs")
                protconf.create_StartingPositions("/home/glenn/Documents/attract/rotation.dat",  receptor.get_filenamePdbReduced(), ligand.get_filenamePdbReduced(),filename_dof)
                benchmark.timer_appendStop("Create_dofs")

        pair = ProteinPair( name_pair, receptor, ligand, filename_dof, receptor.get_pathInput(), receptor.get_pathOutput() , filename_modesJoined = filename_modesJoined,  filename_modesJoined_reference=filename_modesJoinedReference)
        pair.filename_pdb_reference_ligand = filename_ligand_reference
        pairs[name_pair] = pair



    if do_minimization:
        print '**************************************************************'
        print "Run Minimization"
        count_minimization = 0
        for  key, pair in pairs.iteritems():
            count_minimization += 1
            receptor = pair.get_receptor()
            ligand = pair.get_ligand()
            print "{}/{} ".format( count_minimization, len(pairs)),"\t--put in docking queue. protein: ", receptor.get_name(), " and ", ligand.get_name()
            filename_output = os.path.join(receptor.get_pathOutput(), receptor.get_name() + filename_extension_dock)
            if not os.path.isfile(filename_output):
                dock.add_ensembleToQueue( id=key, filename_dofs=pair.filename_dof, filename_output=filename_output,    filename_parameter=parameter_dir,
                                           filename_pdbReceptor=receptor.get_filenamePdbReduced(),               filename_alphabetReceptor=receptor.get_filenameAlphabet(),
                                          filename_gridReceptor=receptor.get_filenameGrid(),                    filename_modesReceptor=receptor.get_filenameModes(),
                                          num_modesReceptor=num_modes,      filename_pdbLigand=ligand.get_filenamePdbReduced(), filename_alphabetLigand=ligand.get_filenameAlphabet(),
                                          filename_gridLigand=ligand.get_filenameGrid(), filename_modesLigand=ligand.get_filenameModes(),
                                          num_modesLigand=num_modes, filename_modesJoined= pair.get_filenameModes(), modeForceFac= 0.5 )
            else:
                print filename_output, " exists already\n"
                finisheditems_docking.put(key)
            time.sleep(0.05)





    if not do_minimization:
        for key, pair in pairs.iteritems():
            finisheditems_docking.put(key)

    dock.stop_threads_if_done()
    if do_scoring:
        # print '**************************************************************'
        # print "Run Scoring"
        # for key, pair in pairs.iteritems():
        count_scoring  = 0
        while not dock.is_done() or not finisheditems_docking.empty() :
            if not finisheditems_docking.empty() :

                pairId = finisheditems_docking.get()
                count_scoring += 1
                pair = pairs[pairId]
                receptor = pair.get_receptor()
                ligand = pair.get_ligand()
                print "{}/{} ".format( count_scoring, len(pairs)),"\t--put in scoring queue. protein: ", receptor.get_name(), " and ", ligand.get_name()
                filename_output = os.path.join(receptor.get_pathOutput(), receptor.get_name() + filename_extension_scoring)

                filename_docking = os.path.join(pair.receptor.get_pathOutput(),
                                                pair.receptor.get_name() + filename_extension_dock)
                if not os.path.isfile(filename_output):
                    score.add_ensembleToQueue(  id = pairId,filename_dofs=filename_docking,filename_output=filename_output,filename_parameter=parameter_dir,
                                              filename_pdbReceptor=receptor.get_filenamePdbReduced(),filename_alphabetReceptor=receptor.get_filenameAlphabet(),
                                              filename_gridReceptor=receptor.get_filenameGrid(),filename_modesReceptor=receptor.get_filenameModes(),
                                              num_modesReceptor=num_modes,filename_pdbLigand=ligand.get_filenamePdbReduced(),filename_alphabetLigand=ligand.get_filenameAlphabet(),
                                              filename_gridLigand=ligand.get_filenameGrid(),filename_modesLigand=ligand.get_filenameModes(),
                                              num_modesLigand=num_modes, filename_modesJoined= pair.get_filenameModes(), radius_cutoff=-1, modeForceFac=0.1 )
                else:
                    print filename_output , " exists already \n"
                    finisheditems_scoring.put(key)
            time.sleep(0.05)

    score.stop_threads_if_done()
    if not do_scoring:
        for key, pair in pairs.iteritems():
            finisheditems_scoring.put(key)
    if do_analyse is True:
        #for count, pair in enumerate( pairs ):
        print '**************************************************************'
        print "Run Analysis"
        count_analysis = 0
        while not score.is_done() or not  finisheditems_scoring.empty():

            if not finisheditems_scoring.empty():
                count_analysis += 1
                #allatom_proteins = protein_ensembles_allatom[count].get_ensemblePDB()
                pairId = finisheditems_scoring.get()
                print "Analyse ", pairId
                pair = pairs[pairId]
                filename_scoring = os.path.join(pair.receptor.get_pathOutput(), pair.receptor.get_name() + filename_extension_scoring)
                filename_docking = os.path.join(pair.receptor.get_pathOutput(), pair.receptor.get_name() + filename_extension_dock)

                print "{}/{} ".format( count_analysis, len(pairs)),"\t--analyse protein: ", pair.receptor.get_name(), " and ", pair.ligand.get_name()
                filename_rmsd = os.path.join(pair.output_folder,pair.receptor.name_protein + "-rmsd.result")
                if not os.path.isfile(filename_rmsd):
                    analyse.run_analysis(pair.output_folder,pair.receptor.name_protein, filename_docking, filename_scoring, os.path.join(pair.input_folder, pair.receptor.filename_allAtom),os.path.join(pair.input_folder,pair.ligand.filename_allAtom),os.path.join(pair.input_folder,pair.filename_pdb_reference_ligand) ,filename_modesReceptor = pair.receptor.get_filenameModes(),
                                     filename_modesLigand=pair.ligand.get_filenameModes(),
                         num_modesReceptor=num_modes, num_modesLigand=num_modes, filename_modesJoined = pair.filename_modesJoined_reference,
                         path_attract=os.environ['ATTRACTDIR'], path_attractTools=os.environ['ATTRACTTOOLS'])
                else:
                    print filename_rmsd, " exists already"


            time.sleep(0.05)


    dock.wait_until_done()
    score.wait_until_done()



    benchmark.timer_addTimer( "Minimization",    benchmark_minimization.sumup_and_getTimer() )
    benchmark.timer_addTimer( "Scoring",         benchmark_scoring.sumup_and_getTimer() )
    benchmark_minimization.save_benchmark( os.path.join(path_folder, name_benchmark + '-MinimizationTime.dat'))
    benchmark.save_benchmark( os.path.join(path_folder, name_benchmark + '-time.dat'))


# run_benchmark( "/home/glenn/Documents/benchmark5_attract/", "-for-docking.pdb",name_benchmark = "benchmark_GPU_0Modes_new", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 0, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=True)

path_test = "/home/glenn/cluster/benchmark_attract_test/1AVX"
path = "/home/glenn/cluster/benchmark5_attract/"
###test
#run_benchmark( path_test, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_0modes_nomode", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 0, use_orig= False
#               , num_threads = 1, do_minimization=True, do_scoring=True)


#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_0modes", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 0, use_orig= False
 #              , num_threads = 1, do_minimization=True, do_scoring=True)

run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_5modes_halfEV", create_grid = True, create_modes = True, create_dofs = True,
               create_reduce = True, num_modes = 5, use_orig= False
               , num_threads = 1, do_minimization=True, do_scoring=True)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_3modes", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 3, use_orig= False
#               , num_threads = 1, do_minimization=True, do_scoring=True)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_1modes", create_grid = True, create_modes = True, create_dofs = True,
 #              create_reduce = True, num_modes = 1, use_orig= False
#               , num_threads = 1, do_minimization=True, do_scoring=True)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_1modes", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 1, use_orig= False, num_threads = 1, do_minimization=True, do_scoring=True)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_3modes", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 3, use_orig= False, num_threads = 1, do_minimization=True, do_scoring=True)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_ORI_scorig_50cut_0modes_2", create_grid = True, create_modes = True, create_dofs = True,
 #              create_reduce = True, num_modes = 0, use_orig= True, num_threads = 1, do_minimization=True, do_scoring=True)


#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_ORI_scorig_50cut_5modes_2", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 5, use_orig= True, num_threads = 1, do_minimization=True, do_scoring=True)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_5modes_2", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 5, use_orig= False, num_threads = 1, do_minimization=True, do_scoring=True)

# run_benchmark( "/home/glenn/Documents/benchmark5_06_04/", "forDock.pdb",name_benchmark = "benchmark_Orig_0Modes_new", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 0, use_orig= True
#                , num_threads = 1, do_minimization=True, do_scoring=True)

#run_benchmark( "/home/glenn/Documents/benchmark5_06_04/", "forDock.pdb",name_benchmark = "benchmark_Orig_5Modes_new", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 5, use_orig= True
#                , num_threads = 1, do_minimization=True, do_scoring=True)




# receptor = protconf.ProteinConfiguration( filename_pdb_protein=filename_receptor, name_protein=name_receptor)
#         receptor.set_path_inputFolder( path_inputFolder=os.path.join(ensemble.get_pathEnsemble(), 'input'))
#         receptor.set_path_outputFolder( path_outputFolder=os.path.join(ensemble.get_pathEnsemble(), name_benchmark))
#         receptor.set_filename( filename_reduce= filename_receptor)
#         receptor.set_chain(chain_receptor)
#         receptor.set_num_modes(num_modes)


