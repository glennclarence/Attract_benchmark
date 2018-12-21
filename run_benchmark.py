import  load_pdbs as lpdb
import ProteinConfiguration as protconf
import os
import run_computation as compute
import run_analysis as analyse
from run_analysis import join_modefiles
import benchmark_timer as benchtime
import time
from Queue import Queue
import compare as compModes
import numpy as np

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
       "compute_attract_orig": "/home/glenn/Downloads/attract_fromHP/bin"}
d_binary={"compute_attract_gpu":"/home/glenn/Documents/Masterarbeit/git/gpuATTRACT_2.0/AttractServer_RELEASE",
       "compute_attract_orig":"/home/glenn/Downloads/attract_fromHP/bin/attract"}

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

        self.filename_modesJoined = filename_modesJoined
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




def run_benchmark( path_folder, filename_scheme, name_benchmark, create_grid = False, create_modes = False, create_dofs = False, create_reduce =False, num_modesRec= 0,num_modesLig = 0, orig_docking= False, orig_scoring = False, num_threads = 7,num_threads_scoring = 1, do_analyse = True, do_scoring = True, do_minimization = True , evfactor = 1.0, rcut = -1, f_dof= None,do_configuration = True,
                   scoring_overwrite=False, analyse_overwrite= False, docking_overwrite = False, analyse_mode = False, useHinsen = False, useAllAtom = False, oModes= False):
    pairs = {}
    finisheditems_scoring = Queue()
    finisheditems_docking = Queue()
    finisheditems_analysis = Queue()
    num_modes = max(num_modesRec, num_modesLig)

    benchmark = benchtime.measure_benchmark()
    benchmark.timer_add("Create_reducepdb")
    benchmark.timer_add("Create_modes")
    benchmark.timer_add("Create_grid")
    benchmark.timer_add("Create_dofs")
    benchmark.timer_add("Analysis")
    benchmark_minimization = benchtime.measure_benchmark()
    benchmark_scoring = benchtime.measure_benchmark()
    benchmark_analysis = benchtime.measure_benchmark()
    #create and start the threads for docking
    dock = compute.Worker( filename_attractBinary = d_binary["compute_attract_orig"] if orig_docking else d_binary["compute_attract_gpu"],
                           do_minimization=True,
                           num_threads=num_threads, args=(benchmark_minimization,finisheditems_docking,), use_OrigAttract=orig_docking)
    dock.start_threads()
    # create and start the threads for scoring
    score = compute.Worker( filename_attractBinary = d_binary["compute_attract_orig"] if orig_scoring else d_binary["compute_attract_gpu"],
                            do_scoring=True, num_threads=num_threads_scoring, args=(benchmark_scoring, finisheditems_scoring,),
                            use_OrigAttract=orig_scoring)
    analysis_threads = compute.Worker(
        filename_attractBinary=d_binary["compute_attract_orig"] if orig_scoring else d_binary["compute_attract_gpu"],
        do_scoring=False, num_threads=num_threads_scoring*2, args=(benchmark_analysis, finisheditems_analysis,),
        use_OrigAttract=False, analyse = True)
    analysis_threads.start_threads()
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
        #print d_config["file_config_receptor_reference"]
        #print name_receptor
        #print os.path.dirname(proteins[name_receptor])
        filename_receptor_reference = os.path.join(os.path.dirname(proteins[name_receptor]), d_config["file_config_receptor_reference"])
        filename_receptor = proteins[ name_receptor]
        receptor = configure_protein(path_inputFolder = path_input,     path_outputFolder = path_output, name_protein = name_receptor,
                                     filename_pdb_protein = filename_receptor,       chain = "A",          num_modes = num_modes, filename_reference=filename_receptor_reference)


        name_ligand = lpdb.get_proteinByindex(proteins,d_config['index_protein'], 'l')
        filename_ligand = proteins[name_ligand]
        filename_ligand_reference = os.path.join( os.path.dirname( proteins[ name_ligand] ), d_config["file_config_ligand_reference"] )
        ligand = configure_protein( path_inputFolder=path_input,    path_outputFolder=path_output,  name_protein=name_ligand,
                                    filename_pdb_protein=filename_ligand,      chain="B",             num_modes=num_modes, filename_reference=filename_ligand_reference)
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
                receptor.set_partner( ligand.get_filenamePdbReduced() )
                if useAllAtom:
                    ligand.set_partner(receptor.get_filenamePdbAllAtom())
                    receptor.set_partner(ligand.get_filenamePdbAllAtom())
                benchmark.timer_start("Create_grid")
                receptor.create_grid(overwrite=False, allAtom=useAllAtom)
                ligand.create_grid(overwrite=False, allAtom = useAllAtom)
                benchmark.timer_appendStop("Create_grid")

            if create_modes and num_modes > 0:
                print "\t\t--Create modefiles"
                benchmark.timer_start("Create_modes")

                try:

                    receptor.create_modes(overwrite=False, heavy= True, allAtom=True)
                    ligand.create_modes(overwrite=False, heavy= True, allAtom=True)
                    if useHinsen:
                        receptor.filename_modes = receptor.name_protein + receptor.ext_hinsen
                        ligand.filename_modes = ligand.name_protein + ligand.ext_hinsen
                    elif oModes:
                        receptor.filename_modes = receptor.name_protein + receptor.ext_oModes
                        ligand.filename_modes = ligand.name_protein + ligand.ext_oModes
                    join_modefiles(receptor.get_filenameModes(), ligand.get_filenameModes(), filename_modesJoined)

                    if useAllAtom:
                        join_modefiles(receptor.get_filenameModesAllAtom(), ligand.get_filenameModesAllAtom(),filename_modesJoined_aa)
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
                            filename_modesJoined= filename_modesJoined )

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
                print pair.filename_modesJoined
                if useAllAtom:
                    dock.add_ensembleToQueue(id=key, filename_dofs=pair.filename_dof, filename_output=filename_output,
                                             filename_parameter=parameter_dir,
                                             filename_pdbReceptor=receptor.get_filenamePdbAllAtom(),
                                             filename_alphabetReceptor=receptor.get_filenameAlphabetAllAtom(),
                                             filename_gridReceptor=receptor.get_filenameGridAllAtom(),
                                             filename_modesReceptor=receptor.get_filenameModesAllAtom(),
                                             num_modesReceptor=num_modesRec,
                                             filename_pdbLigand=ligand.get_filenamePdbAllAtom(),
                                             filename_alphabetLigand=ligand.get_filenameAlphabetAllAtom(),
                                             filename_gridLigand=ligand.get_filenameGridAllAtom(),
                                             filename_modesLigand=ligand.get_filenameModesAllAtom(),
                                             num_modesLigand=num_modesLig,
                                             filename_modesJoined=pair.filename_modesJoined_aa,
                                             modeForceFac=evfactor,
                                             logfile=os.path.join(receptor.get_pathOutput(), "logfile_minimization"))
                else:
                    dock.add_ensembleToQueue( id=key, filename_dofs=pair.filename_dof, filename_output=filename_output,    filename_parameter=parameter_dir,
                                           filename_pdbReceptor=receptor.get_filenamePdbReduced(),               filename_alphabetReceptor=receptor.get_filenameAlphabet(),
                                          filename_gridReceptor=receptor.get_filenameGrid(),                    filename_modesReceptor=receptor.get_filenameModes(),
                                          num_modesReceptor=num_modesRec,      filename_pdbLigand=ligand.get_filenamePdbReduced(), filename_alphabetLigand=ligand.get_filenameAlphabet(),
                                          filename_gridLigand=ligand.get_filenameGrid(), filename_modesLigand=ligand.get_filenameModes(),
                                          num_modesLigand=num_modesLig, filename_modesJoined= pair.filename_modesJoined,
                                          modeForceFac= evfactor , logfile=os.path.join(receptor.get_pathOutput(), "logfile_minimization"))
            else:
                print filename_output, " exists already\n"
                finisheditems_docking.put(key)





    if not do_minimization:
        for key, pair in pairs.iteritems():
            finisheditems_docking.put(key)



    if do_scoring:
        # print '**************************************************************'
        # print "Run Scoring"
        # for key, pair in pairs.iteritems():
        count_scoring  = 0
        #while not dock.is_done() or not finisheditems_docking.empty() :
        #while ( dock.is_comp() or dock.is_empty()) or not finisheditems_docking.empty():
        while count_scoring < len(pairs):

            if not finisheditems_docking.empty() :
                time.sleep(0.2)
                pairId = finisheditems_docking.get()
                count_scoring += 1
                pair = pairs[pairId]
                receptor = pair.receptor
                ligand = pair.ligand
                print "{}/{} ".format( count_scoring, len(pairs)),"\t--put in scoring queue. protein: ", receptor.get_name(), " and ", ligand.get_name()
                filename_output = os.path.join(receptor.get_pathOutput(), receptor.get_name() + d_config['file_compute_scoring'])
                filename_docking = os.path.join(receptor.get_pathOutput(), receptor.get_name() + d_config['file_compute_dock'])
                #print filename_docking
                if not os.path.isfile(filename_output) or scoring_overwrite :
                    print "modes joined ",pair.filename_modesJoined
                    if useAllAtom:
                        score.add_ensembleToQueue(id=pairId, filename_dofs=filename_docking,
                                              filename_output=filename_output, filename_parameter=parameter_dir,
                                              filename_pdbReceptor=receptor.get_filenamePdbAllAtom(),
                                              filename_alphabetReceptor=receptor.get_filenameAlphabetAllAtom(),
                                              filename_gridReceptor=receptor.get_filenameGridAllAtom(),
                                              filename_modesReceptor=receptor.get_filenameModesAllAtom(),
                                              num_modesReceptor=num_modesRec,
                                              filename_pdbLigand=ligand.get_filenamePdbAllAtom(),
                                              filename_alphabetLigand=ligand.get_filenameAlphabetAllAtom(),
                                              filename_gridLigand=ligand.get_filenameGridAllAtom(),
                                              filename_modesLigand=ligand.get_filenameModesAllAtom(),
                                              num_modesLigand=num_modesLig,
                                              filename_modesJoined=pair.filename_modesJoined_aa, radius_cutoff=rcut,
                                              modeForceFac=evfactor,
                                              logfile=os.path.join(receptor.get_pathOutput(), "logfile_scoring"))
                    else:
                        score.add_ensembleToQueue(  id = pairId,filename_dofs=filename_docking,filename_output=filename_output,filename_parameter=parameter_dir,
                                              filename_pdbReceptor=receptor.get_filenamePdbReduced(),filename_alphabetReceptor=receptor.get_filenameAlphabet(),
                                              filename_gridReceptor=receptor.get_filenameGrid(),filename_modesReceptor=receptor.get_filenameModes(),
                                              num_modesReceptor=num_modesRec,filename_pdbLigand=ligand.get_filenamePdbReduced(),filename_alphabetLigand=ligand.get_filenameAlphabet(),
                                              filename_gridLigand=ligand.get_filenameGrid(),filename_modesLigand=ligand.get_filenameModes(),
                                              num_modesLigand=num_modesLig, filename_modesJoined= pair.filename_modesJoined, radius_cutoff=rcut, modeForceFac=evfactor, logfile=os.path.join(receptor.get_pathOutput(), "logfile_scoring") )
                else:
                    print filename_output , " exists already\n"
                    finisheditems_scoring.put(pairId)


    if not do_scoring:
        for key, pair in pairs.iteritems():
            finisheditems_scoring.put(key)

    if do_analyse or analyse_mode is True:
        #for count, pair in enumerate( pairs ):
        print '**************************************************************'
        print "Run Analysis"
        count_analysis = 0
        numAnalysed = 0
        #while not done8 or not  finisheditems_scoring.empty() :
        #while (score.is_comp() or score.is_empty())  or not finisheditems_scoring.empty():
        while numAnalysed < len(protein_ensembles):
            #print score.is_comp(), score.is_empty(),score.is_signal()
            if not finisheditems_scoring.empty():
                count_analysis += 1
                #allatom_proteins = protein_ensembles_allatom[count].get_ensemblePDB()
                pairId = finisheditems_scoring.get()
                numAnalysed +=1
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
                if  do_analyse:


                    print "put task in ana list0"

                    task1 = {'id': pairId,
                           'path_analysis': pair.output_folder, 'name_analysis': pair.receptor.name_protein,
                           'filename_dockResult': filename_docking, 'filename_scoreResult': filename_scoring,
                           'filename_pdbReceptor': fn_rec, 'filename_pdbLigand': fn_lig,
                           'filename_pdbReceptor_aa': fn_rec_aa, 'filename_pdbLigand_aa': fn_lig_aa,
                           'filename_pdbReceptor_heavy': fn_rec_heavy, 'filename_pdbLigand_heavy': fn_lig_heavy,
                            'filename_pdbReceptorRef': fn_rec_ref, 'filename_pdbLigandRef': fn_lig_ref,
                            'filename_pdbReceptorRef_aa': fn_rec_ref_aa, 'filename_pdbLigandRef_aa': fn_lig_ref_aa,
                            'filename_pdbReceptorRef_heavy': fn_rec_ref_heavy,
                            'filename_pdbLigandRef_heavy': fn_lig_ref_heavy,
                            'filename_modesJoined': pair.filename_modesJoined,
                            'filename_modesJoined_aa': pair.filename_modesJoined_aa,
                            'filename_modesJoined_heavy': pair.filename_modesJoined_heavy,
                            'num_modesReceptor': num_modesRec, 'num_modesLigand': num_modesLig, 'path_python': '/usr/bin/python2',
                            'path_attract': '/home/glenn/Downloads/attract_fromHP/bin',
                            'path_attractTools': os.environ['ATTRACTTOOLS'], 'overwrite': analyse_overwrite
                           }
                    #task1 ={'id':pairId}
                    print "put task in ana list"
                    analysis_threads.add_analysiseToQueue(task1)
                	# analyse.run_analysis(pair.output_folder, pair.receptor.name_protein, filename_docking, filename_scoring,
                     #                  filename_pdbReceptor=          fn_rec,         filename_pdbLigand=           fn_lig,
                     #                  filename_pdbReceptor_aa=       fn_rec_aa,      filename_pdbLigand_aa=        fn_lig_aa,
                     #                  filename_pdbReceptor_heavy=    fn_rec_heavy,   filename_pdbLigand_heavy=     fn_lig_heavy,
                     #                  filename_pdbReceptorRef=       fn_rec_ref,     filename_pdbLigandRef=        fn_lig_ref,
                     #                  filename_pdbReceptorRef_aa=    fn_rec_ref_aa,  filename_pdbLigandRef_aa=     fn_lig_ref_aa,
                     #                  filename_pdbReceptorRef_heavy= fn_rec_ref_heavy,   filename_pdbLigandRef_heavy=  fn_lig_ref_heavy,
                     #                  filename_modesJoined=pair.filename_modesJoined, filename_modesJoined_aa=pair.filename_modesJoined_aa,
                     #                  filename_modesJoined_heavy=pair.filename_modesJoined_heavy,num_modesReceptor=num_modesRec, num_modesLigand=num_modesLig,
                     #                  path_attract='/home/glenn/Downloads/attract_fromHP/bin', path_attractTools=os.environ['ATTRACTTOOLS'], overwrite=analyse_overwrite)
                folder_mode = os.path.join(os.path.dirname(pair.output_folder), "mode_eval")

                if analyse_mode:


                    os.system("mkdir -p {} ".format(folder_mode))
                    modes_vec = np.linspace(1, num_modes-1, num_modes -2,dtype=int)#
		    #print modes_vec
                    #print pair.receptor.get_filenameModes()
                    out_pdbs_rec =[]
                    out_pdbs_lig = []
                    for mode in modes_vec:
                        out_pdbs_rec.append(os.path.join(folder_mode, "out_rec_{}.pdb".format(mode)))
                        out_pdbs_lig.append(os.path.join(folder_mode, "out_lig_{}.pdb".format(mode)))

	                compModes.evaluateModes(modes_vec,pair.receptor.get_filenameModes(),
                                        os.path.join(pair.input_folder, pair.receptor.filename_reduce_ref),fn_rec,"dof.dat",  out_pdbs_rec ,os.path.join(folder_mode,"result_rec.dat") )
                    compModes.evaluateModes(modes_vec, pair.ligand.get_filenameModes(),
                                            os.path.join(pair.input_folder, pair.ligand.filename_reduce_ref), fn_lig,
                                            "dof.dat", out_pdbs_lig, os.path.join(folder_mode, "result_lig.dat"))

            else:
                time.sleep(1)

    dock.stop_threads_if_done()
    score.stop_threads_if_done()
    print " finished putting items in queues"
    analysis_threads.stop_threads_if_done()
    print "give stop signal to analyse thread"
    dock.wait_until_done()
    print "stopped docking threads"
    score.wait_until_done()
    print "stopped scoring threads"
    analysis_threads.wait_until_done()
    print "stopped analyse threads"

    benchmark.timer_addTimer( "Minimization",    benchmark_minimization.sumup_and_getTimer() )
    benchmark.timer_addTimer( "Scoring",         benchmark_scoring.sumup_and_getTimer() )
    benchmark.timer_addTimer("Analysing", benchmark_analysis.sumup_and_getTimer())
    benchmark_minimization.save_benchmark( os.path.join(path_folder, name_benchmark + '-MinimizationTime.dat'))
    benchmark_scoring.save_benchmark(os.path.join(path_folder, name_benchmark + '-ScoringTime.dat'))
    benchmark.save_benchmark( os.path.join(path_folder, name_benchmark + '-time.dat'))


# run_benchmark( "/home/glenn/Documents/benchmark5_attract/", "-for-docking.pdb",name_benchmark = "benchmark_GPU_0Modes_new", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modes = 0, use_orig= False
#                , num_threads = 1, do_minimization=True, do_scoring=True)

path_test = "/home/glenn/work/benchmark_attract_test"
path = "/home/glenn/work/benchmark5_attract"
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


#new  docking run 180824
#run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr0_ml0_ev1p0_sO_c50_mr0_ml0_ev1", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modesRec = 0,num_modesLig = 0, orig_docking= False, orig_scoring=True, rcut=50,
#                num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)


# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr5_ml5_ev0p1_sO_c50_mr5_ml5_ev0p1", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 5,num_modesLig = 5, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 0.1, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr5_ml5_ev0p5_sO_c50_mr5_ml5_ev0p5", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 5,num_modesLig = 5, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 0.5, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr5_ml5_ev2p0_sO_c50_mr5_ml5_ev2p0", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 5,num_modesLig = 5, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 2.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr1_ml1_ev1p0_sO_c50_mr1_ml1_ev1p0", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 1,num_modesLig = 1, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)

# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr1_ml1_ev0p5_sO_c50_mr1_ml1_ev0p5", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 1,num_modesLig = 1, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 0.5, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr1_ml1_ev2p0_sO_c50_mr1_ml1_ev2p0", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 1,num_modesLig = 1, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 2.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr1_ml1_ev5p0_sO_c50_mr1_ml1_ev5p0", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 1,num_modesLig = 1, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 5.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr1_ml1_ev0p1_sO_c50_mr1_ml1_ev0p1", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 1,num_modesLig = 1, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 0.1, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr3_ml3_ev1p0_sO_c50_mr3_ml3_ev1p0", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 3,num_modesLig = 3, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr10_ml10_ev1p0_sO_c50_mr10_ml10_ev1p0", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 10,num_modesLig = 10, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr10_ml10_ev2p0_sO_c50_mr10_ml10_ev2p0", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 10,num_modesLig = 10, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 2.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr10_ml10_ev5p0_sO_c50_mr10_ml10_ev5p0", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 10,num_modesLig = 10, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 5.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr10_ml10_ev0p5_sO_c50_mr10_ml10_ev0p5", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 10,num_modesLig = 10, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 0.5, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr10_ml10_ev0p1_sO_c50_mr10_ml10_ev0p1", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 10,num_modesLig = 10, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 0.1, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr15_ml15_ev1p0_sO_c50_mr15_ml15_ev1p0", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 15,num_modesLig = 15, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)

import math

# path = "/home/glenn/work/benchmark5_best"
# numModesList = [15,20]
# numScales = [0.3,1.5,3.0,7.0]
# for modes in numModesList:
#     for scale in numScales:
#         frac, whole = math.modf(scale)
#
#         run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr{}_ml{}_ev{}p{}_sO_c50_mr{}_ml{}_ev{}p{}".format(modes,modes,int(whole),str(frac)[2], modes,modes,int(whole),str(frac)[2]), create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = modes,num_modesLig = modes, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = scale, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
#
# path = "/home/glenn/work/benchmark5_worst"
# numModesList = [1,3,5,10,15,20]
# numScales = [0.3,1.5,3.0,7.0]
# for modes in numModesList:
#     for scale in numScales:
#         frac, whole = math.modf(scale)
#
#         run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr{}_ml{}_ev{}p{}_sO_c50_mr{}_ml{}_ev{}p{}".format(modes,modes,int(whole),str(frac)[2], modes,modes,int(whole),str(frac)[2]), create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = modes,num_modesLig = modes, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = scale, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)


# path = "/home/glenn/Documents/WebNma/best"
# numModesListRec = [1,5,10]
# numModesListLig = [1,5,10]
# numScales = [0.01,0.05,0.1,1]
#
# benchmark_useHinsen = True
# benchmark_allAtom = False
# extention =""
# if benchmark_useHinsen:
#     extention += "_hinsen"
# if benchmark_allAtom:
#     extention += "_aa"
# for modesRec in numModesListRec:
#     for modesLig in numModesListLig:
#         for scale in numScales:
#             if modesLig == modesRec:
#                 frac, whole = math.modf(scale)
#                 run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr{}_ml{}_ev{}p{}_sO_c50_mr{}_ml{}_ev{}p{}{}".format(modesRec,modesLig,int(whole),str(frac)[2:4], modesRec,modesLig,int(whole),str(frac)[2:4],extention), create_grid = True, create_modes = True, create_dofs = True,
#                  create_reduce = True, num_modesRec = modesRec,num_modesLig = modesLig, orig_docking= False, orig_scoring=True, rcut=50,
#                  num_threads = 1, do_minimization=True, do_scoring=True, evfactor = scale, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False, useHinsen=benchmark_useHinsen, useAllAtom=benchmark_allAtom)
#


# pathList =[ "/home/glenn/Documents/WebNma/worst","/home/glenn/Documents/WebNma/best"]
# pathList= ["/home/glenn/work/benchmark5_best","/home/glenn/work/benchmark5_worst"]
# numModesListRec = [0,1,3,5,10,15,20]
# numModesListLig = [0,1,3,5,10,15,20]
#
# logfile = open('logFile_run180924.log', 'w+')
# benchmark_useHinsen = False
# benchmark_allAtom = False
# modeList = [True, False]
# for i in modeList:
#     extention =""
#     benchmark_useHinsen = i
#     if benchmark_useHinsen:
#         extention += "_hinsen"
#         numScales = [0.0001, 0.001, 0.01, 0.05, 0.1, 1]
#     else:
#         numScales = [0.1,0.5,1,2,5]
#     if benchmark_allAtom:
#         extention += "_aa"
#     for path in pathList:
#         for modesRec in numModesListRec:
#             for modesLig in numModesListLig:
#                 for scale in numScales:
#                     if modesRec != 0 or modesLig != 0 :
#                     #if modesLig == modesRec:
#                         frac, whole = math.modf(scale)
#                         try:
#                             run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr{}_ml{}_ev{}p{}_sO_c50_mr{}_ml{}_ev{}p{}{}".format(modesRec,modesLig,int(whole),str(frac)[2:5], modesRec,modesLig,int(whole),str(frac)[2:5],extention), create_grid = True, create_modes = True, create_dofs = True,
#                             create_reduce = True, num_modesRec = modesRec,num_modesLig = modesLig, orig_docking= False, orig_scoring=True, rcut=50,
#                             num_threads = 1, do_minimization=True, do_scoring=True, evfactor = scale, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False, useHinsen=benchmark_useHinsen, useAllAtom=benchmark_allAtom)
#                         except:
#                             logfile.write("failed at  mr {} ml {} scale {} hinsen {}  path {}\n".format( modesRec , modesLig, scale, benchmark_useHinsen, path))
#                             pass

#optimal modes docking
pathList =[ "/home/glenn/Documents/WebNma/worst","/home/glenn/Documents/WebNma/best"]
pathList= ["/home/glenn/work/benchmark5_best","/home/glenn/work/benchmark5_worst"]
pathList= ["/home/glenn/work/benchmark5_attract"]

numModesListRec = [0,1]
numModesListLig = [0,1]
#numModesListLig = [0]
#numModesListRec = [0]
logfile = open('logFile_run181126_2.log', 'w+')
benchmark_useHinsen = True
benchmark_allAtom = False
omodes = False
modeList = [False]
numScales=[   0.00015,0.00018, 0.0002]#0.00008,0.000080.0002,0.0003,0.0004,0.0005,0.0006,0.001,0.002,0.005,0.01,0.02,0.1,1  #0.0015
#numScales=[1.0]
extention = "_hinsen"
for path in pathList:
    for modesRec in numModesListRec:
        for modesLig in numModesListLig:
            for scale in numScales:
                if modesLig == 0 and  modesRec == 0:
                    continue
                frac, whole = math.modf(scale)
                try:
                    run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr{}_ml{}_ev{}p{}_sO_c50_mr{}_ml{}_ev{}p{}{}".format(modesRec,modesLig,int(whole),"{:6f}".format(frac)[2:7], modesRec,modesLig,int(whole),"{:6f}".format(frac)[2:7],extention), create_grid = True, create_modes = True, create_dofs = True,
                    create_reduce = True, num_modesRec = modesRec,num_modesLig = modesLig, orig_docking= False, orig_scoring=True, rcut=50,
                    num_threads = 1, do_minimization=True, do_scoring=True, evfactor = scale, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False, useHinsen=benchmark_useHinsen,oModes=omodes, useAllAtom=benchmark_allAtom,num_threads_scoring = 2)
                except:
                    logfile.write("failed at  mr {} ml {} scale {} hinsen {}  path {}\n".format( modesRec , modesLig, scale, benchmark_useHinsen, path))
                    pass

numModesListRec = [0,3]
numModesListLig = [0,3]
for path in pathList:
    for modesRec in numModesListRec:
        for modesLig in numModesListLig:
            for scale in numScales:
                if modesLig == 0 and  modesRec == 0:
                    continue
                frac, whole = math.modf(scale)
                try:
                    run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr{}_ml{}_ev{}p{}_sO_c50_mr{}_ml{}_ev{}p{}{}".format(modesRec,modesLig,int(whole),"{:6f}".format(frac)[2:7], modesRec,modesLig,int(whole),"{:6f}".format(frac)[2:7],extention), create_grid = True, create_modes = True, create_dofs = True,
                    create_reduce = True, num_modesRec = modesRec,num_modesLig = modesLig, orig_docking= False, orig_scoring=True, rcut=50,
                    num_threads = 1, do_minimization=True, do_scoring=True, evfactor = scale, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False, useHinsen=benchmark_useHinsen,oModes=omodes, useAllAtom=benchmark_allAtom,num_threads_scoring = 2)
                except:
                    logfile.write("failed at  mr {} ml {} scale {} hinsen {}  path {}\n".format( modesRec , modesLig, scale, benchmark_useHinsen, path))
                    pass

numModesListRec = [0,10]
numModesListLig = [0,10]
for path in pathList:
    for modesRec in numModesListRec:
        for modesLig in numModesListLig:
            for scale in numScales:
                if modesLig == 0 and  modesRec == 0:
                    continue
                frac, whole = math.modf(scale)
                try:
                    run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr{}_ml{}_ev{}p{}_sO_c50_mr{}_ml{}_ev{}p{}{}".format(modesRec,modesLig,int(whole),"{:6f}".format(frac)[2:7], modesRec,modesLig,int(whole),"{:6f}".format(frac)[2:7],extention), create_grid = True, create_modes = True, create_dofs = True,
                    create_reduce = True, num_modesRec = modesRec,num_modesLig = modesLig, orig_docking= False, orig_scoring=True, rcut=50,
                    num_threads = 1, do_minimization=True, do_scoring=True, evfactor = scale, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False, useHinsen=benchmark_useHinsen,oModes=omodes, useAllAtom=benchmark_allAtom,num_threads_scoring = 2)
                except:
                    logfile.write("failed at  mr {} ml {} scale {} hinsen {}  path {}\n".format( modesRec , modesLig, scale, benchmark_useHinsen, path))
                    pass


#docking with bound structures
# pathList =[ "/home/glenn/Documents/WebNma/worst","/home/glenn/Documents/WebNma/best"]
# pathList= ["/home/glenn/work/benchmark5_best","/home/glenn/work/benchmark5_worst"]
# pathList= ["/home/glenn/work/benchmark5_bound"]
# numModesListRec = [0]

# numModesListLig = [0]
#
# logfile = open('logFile_run180924.log', 'w+')
# benchmark_useHinsen = False
# benchmark_allAtom = False
# omodes = True
# modeList = [False]
# numScales=[1]
#
# extention = "_bound"
# for path in pathList:
#     for modesRec in numModesListRec:
#         for modesLig in numModesListLig:
#             for scale in numScales:
#
#                 #if modesLig == modesRec:
#                 frac, whole = math.modf(scale)
#                 try:
#                     run_benchmark( path, "-refe.pdb",name_benchmark = "dG_mr{}_ml{}_ev{}p{}_sO_c50_mr{}_ml{}_ev{}p{}{}".format(modesRec,modesLig,int(whole),str(frac)[2:5], modesRec,modesLig,int(whole),str(frac)[2:5],extention), create_grid = True, create_modes = True, create_dofs = True,
#                     create_reduce = True, num_modesRec = modesRec,num_modesLig = modesLig, orig_docking= False, orig_scoring=True, rcut=50,
#                     num_threads = 1, do_minimization=True, do_scoring=True, evfactor = scale, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False, useHinsen=benchmark_useHinsen,oModes=False, useAllAtom=benchmark_allAtom)
#                 except:
#                     logfile.write("failed at  mr {} ml {} scale {} hinsen {}  path {}\n".format( modesRec , modesLig, scale, benchmark_useHinsen, path))
#                     pass
#


#path = "/home/glenn/work/test_independet"
#run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr10_ml10_ev0p3_sO_c50_mr10_ml10_ev0p3_new", create_grid = True, create_modes = True, create_dofs = True,
#                            create_reduce = True, num_modesRec =10,num_modesLig = 10, orig_docking= False, orig_scoring=True, rcut=50,
#                            num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 0.3, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False, useHinsen=False, useAllAtom=False)



#path = "/home/glenn/work/test_hinsen"
#run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr{}_ml{}_ev{}p{}_sO_c50_mr{}_ml{}_ev{}p{}_1".format(5,0,1,0, 5,0,1,0), create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modesRec = 5,num_modesLig =5, orig_docking= False, orig_scoring=True, rcut=50,
#                num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1, do_analyse = True, scoring_overwrite=True, analyse_overwrite=True, docking_overwrite=True,analyse_mode = True, useHinsen=False)


#run_benchmark( path, "-for-docking.pdb",name_be# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr5_ml5_ev0p1_sO_c50_mr5_ml5_ev0p1", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 5,num_modesLig = 5, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 0.1, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr5_ml5_ev0p5_sO_c50_mr5_ml5_ev0p5", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 5,num_modesLig = 5, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 0.5, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr5_ml5_ev2p0_sO_c50_mr5_ml5_ev2p0", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 5,num_modesLig = 5, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 2.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
#
#
# run_benchmark( path, "-for-docking.pdb",name_benchmark = "dG_mr1_ml1_ev1p0_sO_c50_mr1_ml1_ev1p0", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 1,num_modesLig = 1, orig_docking= False, orig_scoring=True, rcut=50,
#                 num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)
# nchmark = "dG_mr20_ml20_ev1p0_sO_c50_mr20_ml20_ev1p0", create_grid = True, create_modes = True, create_dofs = True,
#                 create_reduce = True, num_modesRec = 20,num_modesLig = 20, orig_docking= False, orig_scoring=True, rcut=50,
#                num_threads = 1, do_minimization=True, do_scoring=True, evfactor = 1.0, do_analyse = True, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = False)

#end

#evaluate modes 180824
#path= "/home/glenn/work/benchmark_attract_test"
#run_benchmark( path, "-for-docking.pdb",name_benchmark = "bn_mode_eval", create_grid = True, create_modes = True, create_dofs = True,
#                create_reduce = True, num_modesRec = 20,num_modesLig = 20, orig_docking= False, orig_scoring=True, rcut=50,
#                num_threads = 1, do_minimization=False, do_scoring=False, evfactor = 1.0, do_analyse = False, scoring_overwrite=False, analyse_overwrite=False, docking_overwrite=False,analyse_mode = True)

#end







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
