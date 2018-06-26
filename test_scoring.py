

import pandas as pd
import numpy as np
import os

path_debug = "/home/glenn/Documents/test_attract"
path_output= os.path.join(path_debug,"output")#
path_reference= os.path.join(path_debug,"reference/1AVX")
path_input = "/home/glenn/Documents/test_attract/reference/1AVX_input"


def getDebugFile( use_modes, use_gpu, scoring):
    ext = ""
    if use_gpu:
        ext += "GPU"
    else:
        ext += "CPU"
    if use_modes:
        ext += "_Modes"
    if scoring:
        ext += "/result.dat"
    else:
        ext += "/dof.dat"
    return  os.path.join(path_output,ext)

def getReferenceFile( use_modes,scoring):
    ext = ""
    if use_modes:
        ext += "Modes"
    else:
        ext += "noModes"
    if scoring:
        ext += "/result.dat"
    else:
        ext += "dof.dat"
    return  os.path.join(path_reference,ext)



def loadData( filename):
    return  pd.read_csv( filename ,sep="\s+;\s+|\s;\s|\s+", engine='python', dtype={'a': np.float32} )


def getMean( data):
    mean = []
    for column in data:
        mean.append(data[column].mean())
    return mean

def compare( list1, list2, atol=1e-05, rtol=1e-07):
    idx_line = 0
    count = 0
    std_dev = 0
    for l1, l2 in zip(list1, list2):
        idx_line += 1
        std_dev += (l1- l2)*(l1- l2)
        if np.isclose(l1, l2, rtol=rtol, atol=atol, equal_nan=False) == False:
            count +=0
            #print(idx_line, l1, l2)
        else:
            count += 1
   # print "accuracy", count/idx_line
    return np.sqrt(std_dev/len(list1))



def run_attract( use_orig, filename_dof, folder_input, do_scoring ,num_modes):

    ext_grid ="grid.grid"
    lig ="ligand"
    rec ="receptor"
    ext_modesRec = "modes"+ rec + ".dat"
    ext_modesLig = "modes" + lig + ".dat"

    ext_pdb = "r.pdb"
    ext_alphabet ="grid.alphabet"
    bash_command = ""
    if not use_orig:

        bash_command += "/home/glenn/Documents/Masterarbeit/git/gpuATTRACT_2.0/AttractServer_DEBUG "
        if do_scoring:
            bash_command += "sc "
        else:
            bash_command += "em "
        bash_command += "--dof "
        bash_command += filename_dof
        bash_command+=" -p "
        bash_command+="${ATTRACTDIR}/../attract.par "
        bash_command+="-r "

        bash_command+=os.path.join( folder_input, rec+ext_pdb)
        bash_command+=" -l "
        bash_command+=os.path.join(folder_input, lig + ext_pdb)

        if num_modes > 0:
            bash_command+=" --numModes "
            bash_command+=str(num_modes)
            bash_command+=" --modesr "
            bash_command+=os.path.join(folder_input, ext_modesRec)
            bash_command+=" --modesl "
            bash_command+=os.path.join(folder_input, ext_modesLig)
        bash_command+=" --alphabetrec "
        bash_command+=os.path.join(folder_input, rec+ext_alphabet)
        bash_command+=" --alphabetlig "
        bash_command+=os.path.join(folder_input, lig + ext_alphabet)
        bash_command+=" --gridrec "
        bash_command+=os.path.join(folder_input, rec + ext_grid)
        bash_command+=" --gridlig "
        bash_command+=os.path.join(folder_input, lig + ext_grid)
    else:
        bash_command += "/home/glenn/Documents/attract/bin/attract "
        bash_command += filename_dof
        bash_command += " ${ATTRACTDIR}/../attract.par "
        bash_command += os.path.join(folder_input, rec + ext_pdb) + " "
        bash_command += os.path.join(folder_input, lig + ext_pdb) + " "
        bash_command += " --fix-receptor "
        bash_command += " --grid 1 "
        bash_command += os.path.join(folder_input, rec + ext_grid)
        if num_modes > 0:
            bash_command += " --grid 2 "
            bash_command += os.path.join(folder_input, lig + ext_grid)
            bash_command += " --modes "
            bash_command += folder_input+ "/hm-all.dat "

        bash_command += " --vmax 1000 "
        if do_scoring:
            bash_command += " --score "
            bash_command += " > " + getReferenceFile( num_modes > 0,do_scoring)

    print  bash_command
    os.system(bash_command)

def test_scoring(use_modes, use_gpu):
    file_ref = getReferenceFile(use_modes, use_gpu)
    data_ref = loadData(file_ref)
    print data_ref
    #mean_ref = mean(data_ref)
    #compare(mean_debug, mean_ref, atol=1e-05, rtol=1e-07)

def test_result(use_modes, use_gpu, scoring):
    file_debug = getDebugFile(use_modes, use_gpu, True)
    data_debug = loadData( file_debug )#
    #ean_debug = getMean(data_debug)

    file_ref = getReferenceFile(use_modes, scoring)
    data_ref = loadData(file_ref)
    print type(data_ref)
    template = "{0:12}|{1:3}{2:13.6g}|{3:3}{4:17.6g}|{5:7}{6:10.3g}|{7:7}{8:10.5g}|{9:6}{10:3}"
    for column in data_debug:
        ratio = compare(data_debug[column], data_ref[column], atol=1e-05, rtol=1e-07)
        print template.format( column,"deb", data_debug[column].mean(),"ref",data_ref[column].mean(), "deb-ref ",data_debug[column].mean() - data_ref[column].mean(),"deb/ref",1-data_debug[column].mean()/ data_ref[column].mean(),"ratio", ratio )
        #print "\tmean deb" ,data_debug[column].mean(), "\tref",data_ref[column].mean(), "\trerr out/ref",  data_debug[column].mean()/ data_ref[column].mean(), "\t aerr out-ref ", data_debug[column].mean() - data_ref[column].mean(),"\tratio", ratio


        #print column
    #for column in data_ref:
        # print data_debug[column], data_ref[column]
       # print column
    #mean_ref = getMean(data_ref)

    #compare(mean_debug, mean_ref, atol=1e-05, rtol=1e-07)


#path_input = "/home/glenn/Documents/Masterarbeit/testfolder/1AVX/input"
file_dof ="/dof_docked_all.dat"
#run_attract( True, path_input + file_dof, path_input, do_scoring=True,num_modes=5)
run_attract( False, path_input + file_dof,path_input, do_scoring=True,num_modes=5)
use_modes = True
use_gpu = True
scoring = True
test_result(use_modes, use_gpu, scoring)