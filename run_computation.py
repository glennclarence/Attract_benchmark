import subprocess
from Queue import Queue
from threading import Thread
import os
import shlex

import consumer_producer as conpro

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
    def __init__(self,  path_attract=os.environ['ATTRACTDIR'], name_attractBinary = "AttractServer", num_threads = 1, do_minimization = False, do_scoring = False, use_OrigAttract = False, args = None ):
        self.num_threads = num_threads
        self.do_minimization = do_minimization
        self.do_scoring = do_scoring
        self.queue= Queue()
        self.path_attract = path_attract
        self.filename_attract = os.path.join( path_attract, name_attractBinary )
        self.name_attractBinary = name_attractBinary
        self.consumer = []
        self.use_origAttract = use_OrigAttract

        for i in range( num_threads ):
            self.consumer.append( conpro.ConsumerThread(self.queue, name='Thread ' + str(i), target=self.consumer_compute, args= args))
        self.producer = conpro.ProducerThread( self.queue, name='producer')
        self.producer.start()

    def consumer_compute(self, queue_item, threadId, args = None):
            args_ensemble = queue_item
            string = "{}  \n".format(threadId)
            print string

            run_program( self.filename_attract, args_ensemble, shell=True)


    def start_threads(self):
        for i in range( self.num_threads ):
            self.consumer[i].start()

    def stop_threads(self):
        for i in range( self.num_threads ):
            self.consumer[i].kill()

    def stop_threads_if_done(self):
        for i in range( self.num_threads ):
            self.consumer[i].stop_if_done()

    def add_ensembleToQueue( self, filename_dofs, filename_parameter, filename_pdbReceptor, filename_pdbLigand,
                            filename_gridReceptor, filename_alphabetReceptor, filename_output,
                            filename_gridLigand=None, filename_alphabetLigand=None,  num_modesReceptor=0, num_modesLigand=0,
                            filename_modesReceptor=None, filename_modesLigand=None, filename_modesJoined = None,  ):

        args = list()
        if not self.use_origAttract:
            if self.do_minimization:
                args.append( " em ")
            elif self.do_scoring:
                args.append( " sc ")
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
        else:
            safeAppend( args, filename_dofs )
            safeAppend( args, filename_parameter )
            safeAppend( args, filename_pdbReceptor )
            safeAppend( args, filename_pdbLigand )
            args.append( " --fix-receptor ")
            if filename_modesJoined is not None:
                args.append(" --modes  ")
                safeAppend( args, filename_modesJoined )
            args.append( " --grid 1  ")
            safeAppend(args, filename_gridReceptor)

            if filename_gridLigand is not None:
                args.append( " --grid 2 ")
                safeAppend( args, filename_gridLigand )
            args.append(" --vmax 1000 ")

            args.append( "> "+filename_output )


        self.producer.add( args )


def run_program( filename_binary, arguments, shell = False ):
    args = list()
    args.append( filename_binary )
    for element  in arguments:
        args.append(element)
    print ' '.join(args)
    #shlex.split(' '.join(args))
    os.system( ' '.join(args))
   # process = subprocess.Popen( split(args), stdout=subprocess.PIPE )
    #while True:
        #output = process.stdout.readline()
        #if output == '' and process.poll() is not None:
        #    break
        #if output:
            #print output.strip()



    # def compute(self, queue):
    #     while True:
    #         args_ensemble = queue.get()
    #         if self.do_minimization:
    #             if use_OrigAttract is False:
    #                 args_ensemble.insert( 0, " em ")
    #             run_program( self.filename_attract, args_ensemble, shell=True)
    #         elif self.do_scoring:
    #             if use_OrigAttract is False:
    #                 args_ensemble.insert(0, " sc ")
    #             else:
    #                 args_ensemble.append( " --score")
    #             run_program(self.filename_attract, args_ensemble)
    #         queue.task_done()
    #
    #
    #
    #
    #
    # def compute_serial(self):
    #     count = 0
    #     while not self.queue.empty():
    #         args_ensemble = self.queue.get()
    #         count += 1
    #         print count
    #         if self.do_minimization:
    #             if self.do_minimization:
    #                 args_ensemble.insert(0, " em ")
    #             run_program(self.filename_attract, args_ensemble, shell=True)
    #         elif self.do_scoring:
    #             if self.do_minimization:
    #                 args_ensemble.insert(0, " sc ")
    #             else:
    #                 args_ensemble.append( " --score")
    #             run_program(self.filename_attract, args_ensemble)
    #
    # def start_Worker(self):
    #     for i in range( self.num_threads ):
    #         worker = Thread(target=self.compute, args=(self.queue,))
    #         worker.setDaemon(True)
    #         worker.start()
    #
    # #def add_ensembleToQueue(self, protein_ensemble ):
    #
    #
    # def run(self):
    #
    #     self.queue.join()
