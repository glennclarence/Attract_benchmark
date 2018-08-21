import  load_pdbs as lpdb
import ProteinConfiguration as protconf
import os
import run_computation as compute
import run_analysis as analyse
from run_analysis import join_modefiles
import benchmark_timer as benchtime
import time
from Queue import Queue


#paths
d_config = {'file_compute_dock':'_dock.result',
                'file_compute_scoring':'_scoring.result',
                'folder_input':'input',
                'file_config_dof':'dofs.dat',
                'file_config_receptor_reference':'receptor-refe.pdb',
                'file_config_ligand_reference':'ligand-refe.pdb',
                'index_protein': 5,
                'index_chain': 5
                }

d_paths={"compute_attract_gpu":"/home/glenn/Documents/Masterarbeit/git/gpuATTRACT_2.0",
       "compute_attract_orig": "/home/glenn/Documents/attract/bin"}

d_binary={"compute_attract_gpu":"/home/glenn/Documents/Masterarbeit/git/gpuATTRACT_2.0/AttractServer",
       "compute_attract_orig":"home/glenn/Documents/attract/bin/attract"}

parameter_dir = os.environ['ATTRACTDIR'] + "/../attract.par"


#path_attract_app = "/home/glenn/Documents/Masterarbeit/git/gpuATTRACT_2.0"
#name_attractBinary = "AttractServer"
#file extensions
#filename_extension_dock = "_dock.result"
#filename_extension_scoring = "_scoring.result"


class ProteinPair:
    def __init__(self,id, receptor, ligand, filename_dof, input_folder, output_folder, filename_modesJoined = None, filename_modesJoined_reference=None, filename_pdb_reference_ligand = None,filename_pdb_reference_receptor = None,filename_receptor_heavy = None, filename_ligand_heavy = None ):
        self.receptor = receptor
        self.ligand = ligand
        self.filename_dof =filename_dof
        self.input_folder = input_folder
        self.output_folder = output_folder
        # self.filename_pdb_reference_ligand = filename_pdb_reference_ligand

        # self.filename_modesJoined_reference = filename_modesJoined_reference
        # self.filename_receptor_heavy = filename_receptor_heavy
        # self.filename_ligand_heavy = filename_ligand_heavy

        self.filename_modesJoined = None
        self.filename_modesJoined_aa  = None
        self.filename_modesJoined_heavy  = None
        self.filename_modesJoined_ref  = None
        self.filename_modesJoined_ref_aa = None
        self.filename_modesJoined_ref_heavy = None

        self.id = id
    # def get_receptor(self):
    #     return self.receptor
    # def get_ligand(self):
    #     return self.ligand
    # def get_filenameDOF(self):
    #     return self.filenameDOF
    # def get_filenameModes(self):
    #     return self.filename_modesJoined


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




def run_benchmark( path_folder, filename_scheme, name_benchmark, create_grid = False, create_modes = False, create_dofs = False, create_reduce =False, num_modes= 0, orig_docking= False, orig_scoring = False, num_threads = 7, do_analyse = True, do_scoring = True, do_minimization = True , evfactor = 1.0, rcut = -1, f_dof= None,do_configuration = True,
                   scoring_overwrite=False, analyse_overwrite= False, docking_overwrite = False):
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
    #create and start the threads for docking
    dock = compute.Worker( filename_attractBinary = d_binary["compute_attract_orig"] if orig_docking else d_binary["compute_attract_gpu"],
                           do_minimization=True,
                           num_threads=num_threads, args=(benchmark_minimization,finisheditems_docking,), use_OrigAttract=orig_docking)
    dock.start_threads()
    # create and start the threads for scoring
    score = compute.Worker( filename_attractBinary = d_binary["compute_attract_orig"] if orig_docking else d_binary["compute_attract_gpu"],
                            do_scoring=True, num_threads=num_threads, args=(benchmark_scoring, finisheditems_scoring,),
                            use_OrigAttract=orig_scoring)
    score.start_threads()


    print '**************************************************************'
    print "Load Protein Pdbs"
    #load all proteins inside one folder into
    protein_ensembles =         lpdb.load_fromFolder( path_folder= path_folder, filename_sheme_pdb=filename_scheme)
    print 'detected',len(protein_ensembles) , "protein ensembles"

    print '**************************************************************'
    print "Create Protein Configurations"
    for count, ensemble in enumerate(protein_ensembles):
        #set paths and get proteins from the ensemble list
        proteins = ensemble.get_ensemblePDB()
        #create the input and output path as subfolder
        path_input =    os.path.join( ensemble.get_pathEnsemble(), d_config['folder_input'])
        path_output =   os.path.join( ensemble.get_pathEnsemble(), name_benchmark)
        filename_dof =  os.path.join( os.path.join(ensemble.get_pathEnsemble(),  d_config['folder_input']),  d_config['file_config_dof'])
        #set the pair name as the name of the last folder level ( assuming it is the pdb-id)
        name_pair = os.path.split( ensemble.get_pathEnsemble() )[-1]



        #get the receptor and the ligand via their size and create Proteinconfigurations from it for receptor and ligand
        name_receptor = lpdb.get_proteinByindex(proteins,d_config['index_protein'], 'r')
        filename_receptor_reference = os.path.join(os.path.dirname(proteins[name_receptor]), d_config["file_config_receptor_reference"])
        filename_receptor = proteins[ name_receptor]
        receptor = configure_protein(path_inputFolder = path_input,     path_outputFolder = path_output, name_protein = name_receptor,
                                     filename_pdb_protein = filename_receptor,       chain = "A",          num_modes = num_modes, filename_reference=filename_receptor)


        name_ligand = lpdb.get_proteinByindex(proteins,d_config['index_protein'], 'l')
        filename_ligand = proteins[name_ligand]
        filename_ligand_reference = os.path.join( os.path.dirname( proteins[ name_ligand] ), d_config["file_config_ligand_reference"] )
        ligand = configure_protein( path_inputFolder=path_input,    path_outputFolder=path_output,  name_protein=name_ligand,
                                    filename_pdb_protein=filename_ligand,      chain="B",             num_modes=num_modes, filename_reference=filename_ligand)
        filename_modesJoined = os.path.join(receptor.get_pathInput(), "-allModes-" + str(num_modes) + ".dat")
        filename_modesJoined_aa = os.path.join(receptor.get_pathInput(), "-allModes-" + str(num_modes) + "-aa.dat")
        filename_modesJoined_heavy = os.path.join(receptor.get_pathInput(), "-allModes-" + str(num_modes) + "-heavy.dat")

        filename_modesJoined_ref = os.path.join(receptor.get_pathInput(), "-allModes-ref-" + str(num_modes) + ".dat")
        filename_modesJoined_ref_aa = os.path.join(receptor.get_pathInput(), "-allModes-ref-" + str(num_modes) + "-aa.dat")
        filename_modesJoined_ref_heavy = os.path.join(receptor.get_pathInput(),"-allModes-ref-" + str(num_modes) + "-heavy.dat")
        heavy = True
        allAtom = True
        if do_configuration:
            print "{}/{} ".format( count + 1, len(protein_ensembles)),"\t--Configure protein: ", name_receptor, " and ", name_ligand
            if create_reduce:
                print "\t\t--Create pdbs"
                benchmark.timer_start("Create_reducepdb")
                ligand.createPDB( overwrite=False, allatom=True, heavy=True, reference=True)
                receptor.createPDB( overwrite=False , allatom=True, heavy = True, reference=True)
                benchmark.timer_appendStop("Create_reducepdb")
            if create_grid:
                print "\t\t--Create grids"
                ligand.set_partner( receptor.get_filenamePdbReduced() )
                receptor.set_partner(ligand.get_filenamePdbReduced())
                benchmark.timer_start("Create_grid")
                receptor.create_grid(overwrite=False)
                ligand.create_grid(overwrite=False)
                benchmark.timer_appendStop("Create_grid")

            if create_modes and num_modes > 0:
                print "\t\t--Create modefiles"
                benchmark.timer_start("Create_modes")

                try:

                    receptor.create_modes(overwrite=False, heavy= True, allAtom=True)
                    ligand.create_modes(overwrite=False, heavy= True, allAtom=True)


                    join_modefiles(receptor.get_filenameModes(), ligand.get_filenameModes(), filename_modesJoined)
                    if heavy:
                        join_modefiles(os.path.join(receptor.get_pathInput(), receptor.filename_modes_heavy) ,os.path.join(receptor.get_pathInput(), ligand.filename_modes_heavy),
                                       filename_modesJoined_heavy)
                    if allAtom:

                        join_modefiles(os.path.join(receptor.get_pathInput(), receptor.filename_modes_aa) ,os.path.join(receptor.get_pathInput(), ligand.filename_modes_aa),
                                       filename_modesJoined_aa)

                except:
                    print "could not create Modefiles"
                    pass
                benchmark.timer_appendStop("Create_modes")

            if create_dofs:
                print "\t\t--Create dofs"
                benchmark.timer_start("Create_dofs")
                protconf.create_StartingPositions("/home/glenn/Documents/attract/rotation.dat",  receptor.get_filenamePdbReduced(), ligand.get_filenamePdbReduced(),filename_dof)
                benchmark.timer_appendStop("Create_dofs")

        pair = ProteinPair( name_pair, receptor, ligand, filename_dof, receptor.get_pathInput(), receptor.get_pathOutput() ,  filename_ligand_heavy=os.path.join( ligand.get_pathInput(), ligand.filename_heavy),
                            filename_receptor_heavy=os.path.join( receptor.get_pathInput(), receptor.filename_heavy),
                            filename_modesJoined= filename_modesJoined if create_modes and num_modes > 0   else None)

        pair.filename_modesJoined_aa = filename_modesJoined_aa

        pair.filename_modesJoined_heavy  = filename_modesJoined_heavy
        pair.filename_pdb_reference_ligand = filename_ligand_reference
        pair.filename_pdb_reference_receptor = filename_receptor_reference
        pairs[name_pair] = pair



    if do_minimization:
        print '**************************************************************'
        print "Run Minimization"
        count_minimization = 0
        for  key, pair in pairs.iteritems():
            count_minimization += 1
            receptor = pair.receptor
            ligand = pair.ligand
            print "{}/{} ".format( count_minimization, len(pairs)),"\t--put in docking queue. protein: ", receptor.get_name(), " and ", ligand.get_name()
            filename_output = os.path.join(receptor.get_pathOutput(), receptor.get_name() + d_config['file_compute_dock'])
            if not os.path.isfile(filename_output) or docking_overwrite:

                dock.add_ensembleToQueue( id=key, filename_dofs=pair.filename_dof, filename_output=filename_output,    filename_parameter=parameter_dir,
                                           filename_pdbReceptor=receptor.get_filenamePdbReduced(),               filename_alphabetReceptor=receptor.get_filenameAlphabet(),
                                          filename_gridReceptor=receptor.get_filenameGrid(),                    filename_modesReceptor=receptor.get_filenameModes(),
                                          num_modesReceptor=num_modes,      filename_pdbLigand=ligand.get_filenamePdbReduced(), filename_alphabetLigand=ligand.get_filenameAlphabet(),
                                          filename_gridLigand=ligand.get_filenameGrid(), filename_modesLigand=ligand.get_filenameModes(),
                                          num_modesLigand=num_modes, filename_modesJoined= pair.filename_modesJoined,
                                          modeForceFac= evfactor , logfile=os.path.join(receptor.get_pathOutput(), "logfile_minimization"))
            else:
                print filename_output, " exists already\n"
                finisheditems_docking.put(key)





    if not do_minimization:
        for key, pair in pairs.iteritems():
            finisheditems_docking.put(key)

    dock.stop_threads_if_done()
    dock.wait_until_done()
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
                receptor = pair.receptor
                ligand = pair.ligand
                print "{}/{} ".format( count_scoring, len(pairs)),"\t--put in scoring queue. protein: ", receptor.get_name(), " and ", ligand.get_name()
                filename_output = os.path.join(receptor.get_pathOutput(), receptor.get_name() + d_config['file_compute_scoring'])
                filename_docking = os.path.join(receptor.get_pathOutput(), receptor.get_name() + d_config['file_compute_dock'])
                print filename_docking
                if not os.path.isfile(filename_output) or scoring_overwrite :
                    print filename_docking
                    score.add_ensembleToQueue(  id = pairId,filename_dofs=filename_docking,filename_output=filename_output,filename_parameter=parameter_dir,
                                              filename_pdbReceptor=receptor.get_filenamePdbReduced(),filename_alphabetReceptor=receptor.get_filenameAlphabet(),
                                              filename_gridReceptor=receptor.get_filenameGrid(),filename_modesReceptor=receptor.get_filenameModes(),
                                              num_modesReceptor=num_modes,filename_pdbLigand=ligand.get_filenamePdbReduced(),filename_alphabetLigand=ligand.get_filenameAlphabet(),
                                              filename_gridLigand=ligand.get_filenameGrid(),filename_modesLigand=ligand.get_filenameModes(),
                                              num_modesLigand=num_modes, filename_modesJoined= pair.filename_modesJoined, radius_cutoff=rcut, modeForceFac=evfactor, logfile=os.path.join(receptor.get_pathOutput(), "logfile_scoring") )
                else:
                    print filename_output , " exists already\n"
                    finisheditems_scoring.put(pairId)


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
                filename_scoring = os.path.join(pair.output_folder, pair.receptor.get_name() + d_config['file_compute_scoring'])
                filename_docking = os.path.join(pair.output_folder, pair.receptor.get_name() + d_config['file_compute_dock'])


                fn_rec = os.path.join(pair.input_folder, pair.receptor.filename_reduce)
                fn_rec_aa = os.path.join(pair.input_folder, pair.receptor.filename_allAtom)
                fn_rec_heavy = os.path.join(pair.input_folder, pair.receptor.filename_heavy)
                fn_rec_ref = pair.receptor.filename_pdb_reference
                fn_rec_ref_aa = os.path.join(pair.input_folder,pair.receptor.filename_allAtom_ref)
                fn_rec_ref_heavy = os.path.join(pair.input_folder,pair.receptor.filename_heavy_ref)

                fn_lig = os.path.join(pair.input_folder, pair.ligand.filename_reduce)
                fn_lig_aa = os.path.join(pair.input_folder, pair.ligand.filename_allAtom)
                fn_lig_heavy = os.path.join(pair.input_folder, pair.ligand.filename_heavy)
                fn_lig_ref = pair.ligand.filename_pdb_reference
                fn_lig_ref_aa = os.path.join(pair.input_folder,pair.ligand.filename_allAtom_ref)
                fn_lig_ref_heavy = os.path.join(pair.input_folder,pair.ligand.filename_heavy_ref)



                print "{}/{} ".format( count_analysis, len(pairs)),"\t--analyse protein: ", pair.receptor.get_name(), " and ", pair.ligand.get_name()
                # analyse.run_analysis(pair.output_folder,pair.receptor.name_protein, filename_docking, filename_scoring,
                #                      filename_refLig,
                #      num_modesReceptor=num_modes, num_modesLigand=num_modes, filename_modesJoined = pair.filename_modesJoined,
                #      path_attract=os.environ['ATTRACTDIR'], path_attractTools=os.environ['ATTRACTTOOLS'], overwrite= False)
                print pair.filename_modesJoined
                analyse.run_analysis(pair.output_folder, pair.receptor.name_protein, filename_docking, filename_scoring,
                                     filename_pdbReceptor=          fn_rec,         filename_pdbLigand=           fn_lig,
                                     filename_pdbReceptor_aa=       fn_rec_aa,      filename_pdbLigand_aa=        fn_lig_aa,
                                     filename_pdbReceptor_heavy=    fn_rec_heavy,   filename_pdbLigand_heavy=     fn_lig_heavy,
                                     filename_pdbReceptorRef=       fn_rec_ref,     filename_pdbLigandRef=        fn_lig_ref,
                                     filename_pdbReceptorRef_aa=    fn_rec_ref_aa,  filename_pdbLigandRef_aa=     fn_lig_ref_aa,
                                     filename_pdbReceptorRef_heavy= fn_rec_ref_heavy,   filename_pdbLigandRef_heavy=  fn_lig_ref_heavy,
                                     filename_modesJoined=pair.filename_modesJoined, filename_modesJoined_aa=pair.filename_modesJoined_aa,
                                     filename_modesJoined_heavy=pair.filename_modesJoined_heavy,num_modesReceptor=num_modes, num_modesLigand=num_modes,
                                     path_attract=os.environ['ATTRACTDIR'], path_attractTools=os.environ['ATTRACTTOOLS'], overwrite=analyse_overwrite)




    score.wait_until_done()



    benchmark.timer_addTimer( "Minimization",    benchmark_minimization.sumup_and_getTimer() )
    benchmark.timer_addTimer( "Scoring",         benchmark_scoring.sumup_and_getTimer() )
    benchmark_minimization.save_benchmark( os.path.join(path_folder, name_benchmark + '-MinimizationTime.dat'))
    benchmark_scoring.save_benchmark(os.path.join(path_folder, name_benchmark + '-ScoringTime.dat'))
    benchmark.save_benchmark( os.path.join(path_folder, name_benchmark + '-time.dat'))


# run_benchmark( "/home/glenn/Documents/benchmark5_attract/", "-for-docking.pdb",name_benchmark = "benchmark_GPU_0Modes_new", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 0, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=True)

path_test = "/home/glenn/cluster/benchmark_attract_test"
#path = "/home/glenn/cluster/benchmark5_attract"
#path = "/home/glenn/cluster/benchmark5_attract_ORIG/"
#path = "/home/glenn/Documents/benchmark_1bgx"

###test
#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_ORIG_scorig_0modes", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 0, use_orig= False
 #              , num_threads = 1, do_minimization=True, do_scoring=True)


#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_0modes", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 0, use_orig= False
 #              , num_threads = 1, do_minimization=True, do_scoring=True)

# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_dof_10k_1", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 5, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=True, f_dof="/home/glenn/cluster/benchmark_1bgx/1BGX/input/dofs.dat")






########STILLL to complete
#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_nocut_1modes", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 1, use_orig= False
#               , num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_nocut_1modes_point5EV", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 1, use_orig= False
#               , num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 0.5)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_nocut_1modes_5EV", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 1, use_orig= False
#               , num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 5.0)


#####end to complete




##from 19.7.2018 - 23.07.2018


# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scGPU_nocut_20modes_1EV", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modes = 20, orig_docking= False, orig_scoring=False,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=True, analyse_overwrite=True, docking_overwrite=True)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scGPU_nocut_20modes_point5EV", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modes = 20, orig_docking= False, orig_scoring=False,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 0.5, do_analyse = True, scoring_overwrite=True, analyse_overwrite=True, docking_overwrite=True)
#
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scGPU_nocut_20modes_2EV", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modes = 20, orig_docking= False, orig_scoring=False,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 2.0, do_analyse = True, scoring_overwrite=True, analyse_overwrite=True, docking_overwrite=True)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scGPU_nocut_20modes_4EV", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 20, orig_docking= False, orig_scoring=False,
#                num_threads = 1, do_minimization=False, do_scoring=True, evfactor = 4.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False)


#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scGPU_nocut_10modes_1EV", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 10, orig_docking= False, orig_scoring=False,
#                num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=True, analyse_overwrite=True, docking_overwrite=True)


#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scGPU_nocut_15modes_1EV", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 15, orig_docking= False, orig_scoring=False,
#                num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=True, analyse_overwrite=True, docking_overwrite=True)

##end










# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_0modes", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes =0, orig_docking= False, orig_scoring=False,
#                num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False)


# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_nocut_1modes", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes =1, orig_docking= False, orig_scoring=False,
#                num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_3modes", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes =3, orig_docking= False, orig_scoring=False,
#                num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_5modes", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes =5, orig_docking= False, orig_scoring=False,
#                num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scGPU_nocut_20modes_point5EV", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes =20, orig_docking= False, orig_scoring=False,
#                num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 0.5, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scGPU_nocut_20modes_2EV", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes =20, orig_docking= False, orig_scoring=False,
#                num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 2.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scGPU_nocut_20modes_4EV", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes =20, orig_docking= False, orig_scoring=False,
#                num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 4.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False)
















#run_benchmark( path_test, "-for-docking.pdb",name_benchmark = "benchmark_ORIG_scorig_5modes", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 5, orig_docking= False, orig_scoring=False
#               , num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0)



#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_nocut_10modes_3EV", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 10, use_orig= False
#               , num_threads = 1, do_minimization=True, do_scoring=True,evfactor = 3)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_nocut_10modes_5EV", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 10, use_orig= False
 #              , num_threads = 1, do_minimization=True, do_scoring=True,evfactor = 5)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_1modes", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 1, use_orig= False, num_threads = 1, do_minimization=True, do_scoring=True)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_3modes", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 3, use_orig= False, num_threads = 1, do_minimization=True, do_scoring=True)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_ORIG_scorig_0modes", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 0, use_orig= True, num_threads = 1, do_minimization=True, do_scoring=True)


#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_ORI_scorig_50cut_5modes_2", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 5, use_orig= True, num_threads = 1, do_minimization=True, do_scoring=True)

#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_GPU_scorig_50cut_5modes_2", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 5, use_orig= False, num_threads = 1, do_minimization=True, do_scoring=True)



#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = 'benchmark_GPU_scGPU_50cut_0modes', create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 0, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=True,evfactor = 1.0, rcut = 50)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = 'benchmark_GPU_scGPU_50cut_3modes', create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 3, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=True,evfactor = 1.0, rcut = 50)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = 'benchmark_GPU_scGPU_50cut_5modes', create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 5, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=True,evfactor = 1.0, rcut = 50)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = 'benchmark_GPU_scGPU_50cut_3modes_scnoEV', create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 3, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=True,evfactor = 0.0, rcut = 50)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = 'benchmark_GPU_scGPU_50cut_5modes_scnoEV', create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 5, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=True,evfactor = 0.0, rcut = 50)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = 'benchmark_GPU_scGPU_50cut_3modes_sc2EV', create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 3, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=True,evfactor = 2.0, rcut = 50)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = 'benchmark_GPU_scGPU_50cut_5modes_sc2EV', create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 5, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=True,evfactor = 2.0, rcut = 50)




############################################################### DOF dependence
#run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_dof_5k", create_grid = True, create_modes = True, create_dofs = True,
#               create_reduce = True, num_modes = 5, use_orig= False
#               , num_threads = 1, do_minimization=True, do_scoring=False, do_analyse=False, f_dof="/home/glenn/Documents/benchmark_1bgx/1BGX/input/dofs_5k.dat")
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_dof_40k", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 5, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=False, do_analyse=False, f_dof="/home/glenn/Documents/benchmark_1bgx/1BGX/input/dofs_40k.dat")
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_dof_50k", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 5, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=False, do_analyse=False, f_dof="/home/glenn/Documents/benchmark_1bgx/1BGX/input/dofs_50k.dat")
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_dof_60k", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 5, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=False, do_analyse=False, f_dof="/home/glenn/Documents/benchmark_1bgx/1BGX/input/dofs_60k.dat")
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_dof_70k", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 5, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=False, do_analyse=False, f_dof="/home/glenn/Documents/benchmark_1bgx/1BGX/input/dofs_70k.dat")
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_dof_80k", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 5, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=False, do_analyse=False, f_dof="/home/glenn/Documents/benchmark_1bgx/1BGX/input/dofs_80k.dat")
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_dof_160k", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 5, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=False, do_analyse=False, f_dof="/home/glenn/Documents/benchmark_1bgx/1BGX/input/dofs_160k.dat")
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "benchmark_dof_320k", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 5, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=False, do_analyse=False, f_dof="/home/glenn/Documents/benchmark_1bgx/1BGX/input/dofs_320k.dat")



#leftovers
#filename_receptor_reference = os.path.join( os.path.dirname( proteins[ name_receptor] ), '1AVX_r_b.pdb' )
# name_receptor = lpdb.get_receptorBySize(   proteins )
# name_receptor = lpdb.get_proteinByindex(proteins, 14, 'r')
#filename_receptor_reference = os.path.join(os.path.dirname(proteins[name_receptor]), '1AVX_r_b.pdb')
#chain_receptor = protconf.get_chainfromName(name_receptor, index_chain)
# chain_ligand = protconf.get_chainfromName(name_ligand, index_chain)
##receptor.set_filename(       filename_reduce= filename_receptor )
#
# filename_modesJoinedReference = os.path.join(receptor.get_pathInput(), "allModes-reference-" + str(num_modes) + ".dat")
#
# if not os.path.isfile(filename_modesJoinedReference):
#     protconf.create_modes(os.path.join(receptor.get_pathInput(), receptor.filename_allAtom), receptor.get_pathInput(),
#                           "/receptor-Modes-reference.dat", num_modes)
#     protconf.create_modes(os.path.join(ligand.get_pathInput(), ligand.filename_allAtom), ligand.get_pathInput(),
#                           "/ligand-Modes-reference.dat", num_modes)
#     # protconf.create_modes(os.path.join(receptor.get_pathInput(), filename_receptor_reference),
#     #                      receptor.get_pathInput(), "/receptor-Modes-reference.dat", num_modes)
#     # protconf.create_modes(os.path.join(ligand.get_pathInput(), filename_ligand_reference),
#     #                      ligand.get_pathInput(), "/ligand-Modes-reference.dat", num_modes)
#     join_modefiles(os.path.join(receptor.get_pathInput(), "receptor-Modes-reference.dat"),
#                    os.path.join(ligand.get_pathInput(), "ligand-Modes-reference.dat"), filename_modesJoinedReference)
# filename_modesJoined = os.path.join(receptor.get_pathInput(),
#                                                        "-allModes-" + str(num_modes) + ".dat")
#                    join_modefiles(receptor.get_filenameModes(), ligand.get_filenameModes(), filename_modesJoined)
