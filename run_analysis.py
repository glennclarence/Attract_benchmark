import sys,os

def run_analysis( path_analysis, name_analysis, filename_dockResult, filename_scoreResult,filename_pdbReceptor, filename_pdbLigand,  filename_pdbReference, filename_modesReceptor = None, filename_modesLigand = None,  num_modesReceptor = 0, num_modesLigand = 0, path_python=sys.executable, path_attract= os.environ['ATTRACTDIR'], path_attractTools = os.environ['ATTRACTTOOLS'], filename_modesJoined = None):


    filename_filled =       os.path.join(path_analysis , name_analysis + "-score.dat")
    filename_sorted =       os.path.join(path_analysis ,name_analysis + "-sorted.dat")
    filename_deredundant =  os.path.join(path_analysis ,name_analysis + "-sorted-dr.dat")
    filename_top =          os.path.join(path_analysis ,name_analysis + "-top.dat")
    filename_demode =       os.path.join(path_analysis ,name_analysis + "-top-demode.dat")
    filename_pdbFinal =     os.path.join(path_analysis ,name_analysis + "-result.pdb")
    filename_rmsd = os.path.join(path_analysis, name_analysis + "-rmsd.result")


    analysis_fillEnergies(path_python, path_attractTools, filename_dockResult, filename_scoreResult, filename_filled )
    bash_command = "{} {}/sort.py {} > {}".format(path_python, path_attractTools, filename_filled, filename_sorted)
   # print bash_command
    os.system( bash_command )
    analysis_removeRedundant(path_attract, filename_sorted, num_modesReceptor, num_modesLigand, path_python,
                             path_attractTools, filename_deredundant)

    analysis_getTop(path_attractTools, filename_deredundant, filename_top)
    if num_modesReceptor > 0 or num_modesLigand > 0:
        analysis_removeModes(path_python, path_attractTools, filename_top, filename_demode)
    if num_modesReceptor > 0 or num_modesLigand > 0:
        analysis_collect(path_attract, filename_demode, filename_pdbReceptor, filename_pdbLigand, filename_pdbFinal)
    else:
        analysis_collect(path_attract, filename_top, filename_pdbReceptor, filename_pdbLigand, filename_pdbFinal)

    calculate_rmsd(path_attract, filename_deredundant, filename_pdbLigand, filename_pdbReference, filename_modesJoined,
                   filename_pdbReceptor,  filename_rmsd, num_modesLigand, num_modesReceptor)
    #print "done1"


def analysis_fillEnergies( path_python, path_attractTools, filename_dockResult, filename_scoreResult, filename_output):
    bash_command = "{} {}/fill-energies.py {} {} > {}".format(path_python, path_attractTools, filename_dockResult,
                                                              filename_scoreResult, filename_output)
   # print bash_command
    os.system(bash_command)

def analysis_removeRedundant( path_attract, filename_sorted, num_modesReceptor, num_modesLigand, path_python, path_attractTools, filename_deredundant):
    bash_command = "{}/deredundant {} 2 --modes {} {} | {} {}/fill-deredundant.py /dev/stdin {} > {}".format(
        path_attract, filename_sorted, num_modesReceptor, num_modesLigand, path_python, path_attractTools,
        filename_sorted, filename_deredundant)
   # print bash_command
    os.system(bash_command)


def analysis_getTop( path_attractTools, filename_deredundant, filename_top, number_best = 50):
    bash_command = "{}/top {} {} > {}".format(path_attractTools, filename_deredundant, number_best, filename_top)
   # print bash_command
    os.system(bash_command)


def analysis_removeModes( path_python, path_attractTools, filename_top, filename_demode):
    bash_command = "{} {}/demode.py {} > {}".format(path_python, path_attractTools, filename_top, filename_demode)
   # print bash_command
    os.system(bash_command)

def analysis_collect( path_attract, filename_demode, filename_pdbReceptor, filename_pdbLigand , filename_pdbFinal):
    bash_command = "{}/collect {} {} {} > {}".format(path_attract, filename_demode, filename_pdbReceptor,
                                                     filename_pdbLigand, filename_pdbFinal)
  #  print bash_command
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
    #print bash_command
    os.system(bash_command)


# bash_command = "{} {}/sort.py {} > {}".format(path_python, path_attractTools, filename_filled, filename_sorted)


  #  bash_command = "{}/deredundant {} 2 -modes {} {} | {} {}/fill-deredundant.py /dev/stdin {} > {}".format(
        #            path_attract, filename_sorted, num_modesReceptor, num_modesLigand, path_python, path_attractTools, filename_sorted, filename_deredundant)


  #  bash_command = "{}/top {} 50 > {}".format( path_attractTools, filename_deredundant, filename_top)

   # bash_command = "{} {}/demode.py {} > {}".format( path_python, path_attractTools, filename_top, filename_demode)

  #  bash_command = "{}/collect {} {} {} > {}".format(path_attract, filename_demode, filename_pdbReceptor, filename_pdbLigand, filename_pdbFinal)