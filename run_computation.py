import subprocess
from Queue import Queue
from threading import Thread
import os
import shlex
#only append a filename to a list if it not None and if it exists
def safeAppend( list, filename):
    if filename is not None and os.path.isfile(filename):
        list.append( filename )
        return True
    else:
        return False
def filecheck(  filename):
    #print filename
    if filename is not None and os.path.isfile(filename):
        return True
    else:
        return False


class Worker:
    def __init__(self,  path_attract=os.environ['ATTRACTDIR'], name_attractBinary = "AttractServer", num_threads = 1, do_minimization = False, do_scoring = False ):
        self.num_threads = num_threads
        self.do_minimization = do_minimization
        self.do_scoring = do_scoring
        self.queue= Queue()
        self.path_attract = path_attract
        self.filename_attract = os.path.join( path_attract, name_attractBinary )
        self.name_attractBinary = name_attractBinary

    def compute(self, queue):
        while True:
            args_ensemble = queue.get()
            if self.do_minimization:
                args_ensemble.insert( 0, " em ")
                run_program( self.filename_attract, args_ensemble, shell=True)
            elif self.do_scoring:
                args_ensemble.insert(0, " sc ")
                run_program(self.filename_attract, args_ensemble)
            queue.task_done()

    def compute_serial(self):
        while not self.queue.empty():
            args_ensemble = self.queue.get()
            #print args_ensemble
            if self.do_minimization:
                args_ensemble.insert(0, " em ")
                run_program(self.filename_attract, args_ensemble, shell=True)
            elif self.do_scoring:
                args_ensemble.insert(0, " sc ")
                run_program(self.filename_attract, args_ensemble)

    def start_Worker(self):
        for i in range( self.num_threads ):
            worker = Thread(target=self.compute, args=(self.queue,))
            worker.setDaemon(True)
            worker.start()

    #def add_ensembleToQueue(self, protein_ensemble ):


    def run(self):

        self.queue.join()

    def add_ensembleToQueue( self, filename_dofs, filename_parameter, filename_pdbReceptor, filename_pdbLigand,
                            filename_gridReceptor, filename_alphabetReceptor, filename_output,
                            filename_gridLigand=None, filename_alphabetLigand=None,  num_modesReceptor=0, num_modesLigand=0,
                            filename_modesReceptor=None, filename_modesLigand=None ):
        args = list()
        args.append(" --dof "), safeAppend(args, filename_dofs)
        args.append(" -p ")
        safeAppend(args, filename_parameter)
        args.append(" --alphabetrec ")
        safeAppend(args, filename_alphabetReceptor)
        if filecheck(filename_alphabetLigand):
            args.append(" --alphabetlig ")
            args.append(filename_alphabetLigand)
        args.append(" --gridrec ")
        safeAppend(args, filename_gridReceptor)
        if filecheck(filename_gridLigand):
            args.append(" --gridlig ")
            safeAppend(args, filename_gridLigand)

        args.append(" -r ")
        safeAppend(args, filename_pdbReceptor)
        args.append(" -l ")
        safeAppend(args, filename_pdbLigand)
        if num_modesReceptor > 0 and filecheck(filename_modesReceptor):
            args.append(" --numModesRec ")
            args.append(str(num_modesReceptor))
            args.append(" --modesr ")
            safeAppend(args, filename_modesReceptor)
        if num_modesLigand > 0 and filecheck(filename_modesLigand):
            args.append(" --numModesLig ")
            args.append(str(num_modesLigand))
            args.append(" --modesl ")
            safeAppend(args, filename_modesLigand)
        args.append(" --output ")
        args.append( filename_output )
        self.queue.put( args )

    #def get_queue(self):
    #    return self.queue

def run_program( filename_binary, arguments, shell = False ):
    args = list()
    args.append( filename_binary )
    for element  in arguments:
        args.append(element)
    #print ''.join(args)
    #os.chdir()
    process = subprocess.Popen( shlex.split(''.join(args)), stdout=subprocess.PIPE )
    #while True:
        #output = process.stdout.readline()
        #if output == '' and process.poll() is not None:
        #    break
        #if output:
            #print output.strip()

    process.wait()

#test_worker = Worker(path_attract="/home/glenn/Documents/Masterthesis/gpuATTRACT_2.0", do_minimization= True)
#test_worker.start_Worker()
#test_worker.run()



#parameter_dir= os.environ['ATTRACTDIR']+"/../attract.par"
#test_worker.add_ensembleToQueue("/home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/dofs.dat",
 #                               parameter_dir,
  #                              "/home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_reduce.pdb",
   #                             "/home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_reduce.pdb",
    #                            "/home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptorgrid.grid",
     #                           "/home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_grid.alphabet",
      #                          "/home/glenn/Documents/Masterthesis/testfolder/benchmark_test/output/output",
       #                         "/home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligandgrid.grid",
        #                        "/home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_grid.alphabet",

#num_modesReceptor = 5, num_modesLigand = 5, filename_modesReceptor="/home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_modes.dat",
#filename_modesLigand ="/home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_modes.dat")


#"--solver"," VA13 ",

#/home/glenn/Documents/Masterthesis/gpuATTRACT_2.0/AttractServer em
# --dof /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/dofs.dat
# -p /home/glenn/Documents/Masterthesis/attract_unchanged/attract.par
# --alphabetrec /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_grid.alphabet
# --alphabetlig /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_grid.alphabet
# --gridrec /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptorgrid.grid
# --gridlig /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligandgrid.grid
# -r /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_reduce.pdb
# -l /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_reduce.pdb
# --numModesRec 5 --modesr /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_modes.dat
# --numModesLig 5 --modesl /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_modes.dat
# > /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/output/output

#/home/glenn/Documents/Masterthesis/gpuATTRACT_2.0/AttractServer em  --dof /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/dofs.dat -p /home/glenn/Documents/Masterthesis/attract_unchanged/attract.par --alphabetrec /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_grid.alphabet --alphabetlig /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_grid.alphabet --gridrec /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptorgrid.grid --gridlig /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligandgrid.grid -r /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_reduce.pdb -l /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_reduce.pdb --numModesRec 5 --modesr /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_modes.dat --numModesLig 5 --modesl /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_modes.dat > /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/output/output
#p=subprocess.call("/home/glenn/Documents/Masterthesis/gpuATTRACT_2.0/AttractServer em  --dof /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/dofs.dat -p /home/glenn/Documents/Masterthesis/attract_unchanged/attract.par --alphabetrec /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_grid.alphabet --alphabetlig /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_grid.alphabet --gridrec /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptorgrid.grid --gridlig /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligandgrid.grid -r /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_reduce.pdb -l /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_reduce.pdb --numModesRec 0 --modesr /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_modes.dat --numModesLig 0 --modesl /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_modes.dat > /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/output/output", shell= True)
#p=subprocess.call("/home/glenn/Documents/Masterthesis/gpuATTRACT_2.0/AttractServer sc  --dof /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/dofs.dat -p /home/glenn/Documents/Masterthesis/attract_unchanged/bin/../attract.par --alphabetrec /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_grid.alphabet  --gridrec /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptorgrid.grid  -r /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/receptor_reduce.pdb -l /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/ligand_reduce.pdb  > /home/glenn/Documents/Masterthesis/testfolder/benchmark_test/output/output", shell= True)