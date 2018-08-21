import sys,os
DEBUG_COMMAND = True
#from ProteinConfiguration import createPDB
#from ProteinConfiguration import create_modes

def run_analysis( path_analysis, name_analysis, filename_dockResult, filename_scoreResult,
                  filename_pdbReceptor, filename_pdbLigand,
                  filename_pdbReceptor_aa, filename_pdbLigand_aa,
filename_pdbReceptor_heavy, filename_pdbLigand_heavy,
                  filename_pdbReceptorRef, filename_pdbLigandRef,
                  filename_pdbReceptorRef_aa, filename_pdbLigandRef_aa,
        filename_pdbReceptorRef_heavy, filename_pdbLigandRef_heavy,
filename_modesJoined = None,filename_modesJoined_aa = None,filename_modesJoined_heavy = None, num_modesReceptor = 0, num_modesLigand = 0, path_python=sys.executable, path_attract= os.environ['ATTRACTDIR'],
                  path_attractTools = os.environ['ATTRACTTOOLS'],  overwrite = False):
    filename_irmsd =        os.path.join(path_analysis, name_analysis + "-irmsd.dat")
    filename_fnat =         os.path.join(path_analysis, name_analysis + "-fnat.dat")
    filename_filled =       os.path.join(path_analysis , name_analysis + "-score.dat")
    filename_sorted =       os.path.join(path_analysis ,name_analysis + "-sorted.dat")
    filename_deredundant =  os.path.join(path_analysis ,name_analysis + "-sorted-dr.dat")
    filename_top =          os.path.join(path_analysis ,name_analysis + "-top.dat")
    filename_demode =       os.path.join(path_analysis ,name_analysis + "-top-demode.dat")
    filename_demode_all = os.path.join(path_analysis, name_analysis + "-top-demode-all.dat")
    filename_pdbFinal =     os.path.join(path_analysis ,name_analysis + "-result.pdb")
    filename_rmsd = os.path.join(path_analysis, name_analysis + "-rmsd.result")

    if not os.path.isfile(filename_filled) or overwrite is True:
        analysis_fillEnergies(path_python, path_attractTools, filename_dockResult, filename_scoreResult, filename_filled )
    if not os.path.isfile(filename_sorted) or overwrite is True:
        bash_command = "{} {}/sort.py {} > {}".format(path_python, path_attractTools, filename_filled, filename_sorted)
        os.system( bash_command )
    if not os.path.isfile(filename_deredundant) or overwrite is True:
        analysis_removeRedundant(path_attract, filename_sorted, num_modesReceptor, num_modesLigand, path_python,
                             path_attractTools, filename_deredundant)
    if not os.path.isfile(filename_top) or overwrite is True:
        analysis_getTop(path_attractTools, filename_deredundant, filename_top)

    if not os.path.isfile(filename_demode) or overwrite is True:
        if num_modesReceptor > 0 or num_modesLigand > 0:
            analysis_removeModes(path_python, path_attractTools, filename_top, filename_demode)


    if num_modesReceptor == 0 and num_modesLigand == 0:
        filename_demode_all = filename_deredundant
        filename_demode = filename_top
    if not os.path.isfile(filename_demode_all) or overwrite is True:
        if num_modesReceptor > 0 or num_modesLigand > 0:
            analysis_removeModes(path_python, path_attractTools, filename_deredundant, filename_demode_all)

    if not os.path.isfile(filename_pdbFinal) or overwrite is True:
        if num_modesReceptor > 0 or num_modesLigand > 0:
            analysis_collect(path_attract, filename_demode, filename_pdbReceptor, filename_pdbLigand, filename_pdbFinal)
        else:
            analysis_collect(path_attract, filename_demode,    filename_pdbReceptor, filename_pdbLigand, filename_pdbFinal)

    if not os.path.isfile( filename_rmsd) or overwrite is True:
        calculate_rmsd(path_attract, filename_demode_all, filename_pdbLigand_aa, filename_pdbLigandRef_aa, filename_modesJoined_aa,
                       filename_pdbReceptor_aa,  filename_rmsd, 0, 0)

    if not os.path.isfile(filename_fnat) or overwrite is True:
        calc_fnat(path_attract, filename_demode,filename_pdbReceptor_heavy,filename_pdbReceptorRef_heavy,
                  filename_pdbLigand_heavy,filename_pdbLigandRef_heavy,filename_fnat, 0, filename_modesJoined_heavy)
    print filename_irmsd, os.path.isfile(filename_irmsd)
    if not os.path.isfile(filename_irmsd) or overwrite is True:
        calc_IRmsd(path_attract, filename_demode, filename_pdbReceptor_heavy, filename_pdbReceptorRef_heavy,
                   filename_pdbLigand_heavy, filename_pdbLigandRef_heavy,filename_irmsd,0,filename_modesJoined_heavy)




def analysis_fillEnergies( path_python, path_attractTools, filename_dockResult, filename_scoreResult, filename_output):
    bash_command = "{} {}/fill-energies.py {} {} > {}".format(path_python, path_attractTools, filename_dockResult,
                                                              filename_scoreResult, filename_output)
   # print bash_command
    if DEBUG_COMMAND:
        print bash_command
    os.system(bash_command)

def analysis_removeRedundant( path_attract, filename_sorted, num_modesReceptor, num_modesLigand, path_python, path_attractTools, filename_deredundant):
    bash_command = "{}/deredundant {} 2 --modes {} {} | {} {}/fill-deredundant.py /dev/stdin {} > {}".format(
        path_attract, filename_sorted, num_modesReceptor, num_modesLigand, path_python, path_attractTools,
        filename_sorted, filename_deredundant)
    if DEBUG_COMMAND:
        print bash_command
    os.system(bash_command)


def analysis_getTop( path_attractTools, filename_deredundant, filename_top, number_best = 50):
    bash_command = "{}/top {} {} > {}".format(path_attractTools, filename_deredundant, number_best, filename_top)
    if DEBUG_COMMAND:
        print bash_command
    os.system(bash_command)


def analysis_removeModes( path_python, path_attractTools, filename_top, filename_demode):
    bash_command = "{} {}/demode.py {} > {}".format(path_python, path_attractTools, filename_top, filename_demode)
    if DEBUG_COMMAND:
        print bash_command
    os.system(bash_command)

def analysis_collect( path_attract, filename_demode, filename_pdbReceptor, filename_pdbLigand , filename_pdbFinal):
    bash_command = "{}/collect {} {} {} > {}".format(path_attract, filename_demode, filename_pdbReceptor,
                                                     filename_pdbLigand, filename_pdbFinal)
    if DEBUG_COMMAND:
        print bash_command
    os.system( bash_command )

def join_modefiles( filename_modesReceptor, filename_modesLigand, filename_output):
    bash_command = "cat /dev/null > {}".format( filename_output)
    os.system( bash_command )
    bash_command = "cat {}  >> {}".format(filename_modesReceptor, filename_output)
    os.system(bash_command)
    bash_command = "cat {}  >> {}".format(filename_modesLigand, filename_output)
    os.system(bash_command)

def calculate_rmsd(path_attract, filename_deredundant , filename_pdbLigandAllAtom, filename_pdbLigandrmsd, modefile,filename_pdbRceptorAllAtom, rmsd_output, num_modesLigand, num_modesReceptor):
    if num_modesReceptor > 0 or num_modesLigand > 0:
        bash_command = "python2 {}/lrmsd.py {} {} {} --modes  {} --receptor {} > {}".format(path_attract, filename_deredundant , filename_pdbLigandAllAtom, filename_pdbLigandrmsd, modefile,filename_pdbRceptorAllAtom
                                                        , rmsd_output)
    else:
        
        bash_command = "python2 {}/lrmsd.py {} {} {}  --receptor {} > {}".format(path_attract,
                                                                                            filename_deredundant,
                                                                                            filename_pdbLigandAllAtom,
                                                                                            filename_pdbLigandrmsd,

                                                                                            filename_pdbRceptorAllAtom
                                                                                            , rmsd_output)
    if DEBUG_COMMAND:
        print bash_command
    #print bash_command
    os.system(bash_command)



def calc_IRmsd(path_attract, filename_demode , filename_pdbReceptorHeavy, filename_pdbReceptorRefe, filename_pdbLigandHeavy, filename_pdbLigandRefe,irmsd_output, num_modes=None, modefile=None):

   # if num_modesReceptor > 0 or num_modesLigand > 0:
    if num_modes > 0:
        bash_command = "python2 {}/irmsd.py {} {} {} {} {} --modes {} > {}".format(path_attract, filename_demode , filename_pdbReceptorHeavy, filename_pdbReceptorRefe,filename_pdbLigandHeavy, filename_pdbLigandRefe,modefile,
                                                         irmsd_output)
    else:
        bash_command = "python2 {}/irmsd.py {} {} {} {} {} > {}".format(path_attract, filename_demode,
                                                                                   filename_pdbReceptorHeavy,
                                                                                   filename_pdbReceptorRefe,
                                                                                   filename_pdbLigandHeavy,
                                                                                   filename_pdbLigandRefe,
                                                                                   irmsd_output)


    if DEBUG_COMMAND:
        print bash_command
    os.system(bash_command)

def calc_fnat(path_attract, filename_demode , filename_pdbReceptorHeavy, filename_pdbReceptorRefe, filename_pdbLigandHeavy, filename_pdbLigandRefe,fnat_output, num_modes=None, modefile=None):
    if num_modes > 0 :
        bash_command = "python2 {}/fnat.py {} 5 {} {} {} {} --modes {} > {}".format(path_attract, filename_demode,
                                                                    filename_pdbReceptorHeavy, filename_pdbReceptorRefe,
                                                                    filename_pdbLigandHeavy, filename_pdbLigandRefe,modefile,
                                                                    fnat_output)
    else:
        bash_command = "python2 {}/fnat.py {} 5 {} {} {} {} > {}".format(path_attract, filename_demode,
                                                                                 filename_pdbReceptorHeavy,
                                                                                 filename_pdbReceptorRefe,
                                                                                 filename_pdbLigandHeavy,
                                                                                 filename_pdbLigandRefe,
                                                                                 fnat_output)
    if DEBUG_COMMAND:
        print bash_command
    os.system(bash_command)


#def run_iAttract():

# bash_command = "{} {}/sort.py {} > {}".format(path_python, path_attractTools, filename_filled, filename_sorted)


  #  bash_command = "{}/deredundant {} 2 -modes {} {} | {} {}/fill-deredundant.py /dev/stdin {} > {}".format(
        #            path_attract, filename_sorted, num_modesReceptor, num_modesLigand, path_python, path_attractTools, filename_sorted, filename_deredundant)


  #  bash_command = "{}/top {} 50 > {}".format( path_attractTools, filename_deredundant, filename_top)

   # bash_command = "{} {}/demode.py {} > {}".format( path_python, path_attractTools, filename_top, filename_demode)

  #  bash_command = "{}/collect {} {} {} > {}".format(path_attract, filename_demode, filename_pdbReceptor, filename_pdbLigand, filename_pdbFinal)