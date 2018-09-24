import subprocess
from Queue import Queue
from threading import Thread
import os
import shlex
import time
import random
from subprocess import PIPE

import consumer_producer as conpro

#only append a filename to a list if it not None and if it exists
def safeAppend( list, filename):
    if filename is not None and os.path.isfile(filename):
        list.append( filename )
        return True
    else:
        return False
def filecheck(  filename):
    if filename is not None and os.path.isfile(filename):
        return True
    else:
        return False


class Worker:
    def __init__(self,  filename_attractBinary = "AttractServer", num_threads = 1, do_minimization = False, do_scoring = False, use_OrigAttract = False, args = None ):
        self.num_threads = num_threads
        self.do_minimization = do_minimization
        self.do_scoring = do_scoring
        self.queue= Queue()
        self.filename_attract = filename_attractBinary
        self.consumer = []
        self.use_origAttract = use_OrigAttract

        for i in range( num_threads ):
            self.consumer.append( conpro.ConsumerThread(self.queue, name='Thread ' + str(i), target=self.consumer_compute, args= args))
        self.producer = conpro.ProducerThread( self.queue, name='producer')
        self.producer.start()

    def consumer_compute(self, queue_item, threadId, args = None):
            task_id = queue_item[0]
            task = queue_item[1]
            #string = "{}  \n".format(threadId)
            #print string
            name_timer = ""
            if self.do_scoring:
                name_timer += "Scoring - "
            elif self.do_minimization:
                name_timer += "Docking - "
            name_timer += str( task_id )
            print  "Started", name_timer
            args[0].timer_add( name_timer, start= True)
            time.sleep(random.random())
            run_program( self.filename_attract, task, shell=True)
            args[0].timer_stop( name_timer )
            print "Finished", name_timer, ". Time: " , args[0].get_elapsedTime(name_timer )
            args[1].put(task_id)

    def get_sizeTask(self):
        return self.queue.qsize()

    def start_threads(self):
        for i in range( self.num_threads ):
            self.consumer[i].start()

    def stop_threads(self):
        for i in range( self.num_threads ):
            self.consumer[i].kill()

    def stop_threads_if_done(self):
        for i in range( self.num_threads ):
            self.consumer[i].stop_if_done()

    def wait_until_done(self):
        done = False
        while not done:
            done = True
            for i in range(self.num_threads):
                if not self.consumer[i].is_done():
                    done = False
            time.sleep( 1 )

    def is_done(self):
        done = True
        for i in range(self.num_threads):
            if not self.consumer[i].is_done():
                done = False
        return done

    def add_ensembleToQueue( self, id, filename_dofs, filename_parameter, filename_pdbReceptor, filename_pdbLigand,
                            filename_gridReceptor, filename_alphabetReceptor, filename_output,
                            filename_gridLigand=None, filename_alphabetLigand=None,  num_modesReceptor=0, num_modesLigand=0,
                            filename_modesReceptor=None, filename_modesLigand=None, filename_modesJoined = None, radius_cutoff = 0, modeForceFac = 1.0 , logfile=None ):

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
            #if (num_modesReceptor > 0 ) or num_modesLigand > 0and filecheck(filename_modesReceptor):
               # args.append(" --numModesRec ")
            args.append(" --numModesRec ")
            args.append(str(num_modesReceptor))
            args.append(" --numModesLig ")
            args.append(str(num_modesLigand))
            args.append(" --modesr ")
            safeAppend(args, filename_modesReceptor)
            if num_modesLigand > 0 and filecheck(filename_modesLigand):
                #args.append(" --numModesLig ")

                #args.append(str(num_modesLigand))
                args.append(" --modesl ")
                safeAppend(args, filename_modesLigand)
            if radius_cutoff > 0 and self.do_scoring:
                args.append("--cutoff ")
                args.append(str( radius_cutoff ))
            if modeForceFac != 1.0:
                args.append( " --evscale  ")
                args.append( str(modeForceFac))

            #args.append("--numCPUs ")
            #args.append(" 16 ")
            #args.append(" --output ")
            args.append(" > ")
            args.append( filename_output )
        else:
            safeAppend( args, filename_dofs )
            safeAppend( args, filename_parameter )
            #print "this is in safe append runcomputatuon filename",filename_pdbReceptor
            safeAppend( args, filename_pdbReceptor )
            safeAppend( args, filename_pdbLigand )
            args.append( " --fix-receptor ")
            if filename_modesJoined is not None and (num_modesLigand > 0 or num_modesReceptor > 0):
                args.append(" --modes  ")
                safeAppend( args, filename_modesJoined )
            #if num_modesReceptor > 0:
            args.append(" --numModesRec  ")
            args.append(str(num_modesReceptor))
            #if num_modesLigand > 0:
            args.append(" --numModesLig  ")
            args.append(str(num_modesLigand))

            #args.append( " --grid 1  ")
            #safeAppend(args, filename_gridReceptor)

            #if filename_gridLigand is not None:
            #    args.append( " --grid 2 ")
            #    safeAppend( args, filename_gridLigand )
            args.append(" --vmax 1000 ")
            if radius_cutoff > 0:
                args.append(" --rcut ")
                args.append(str(radius_cutoff))
            if modeForceFac != 1.0:
                args.append( " --evscale  ")
                args.append( str(modeForceFac))
            if self.do_scoring:
                args.append( " --score ")
            args.append( "> "+filename_output )

        if logfile is not None:
            with open(logfile, 'w') as f:
                f.write(''.join(args))
                f.close()
        self.producer.add( (id,args) )


def run_program( filename_binary, arguments, shell = False ):
    args = list()
    args.append( filename_binary )
    for element  in arguments:
        args.append(element)
    print ' '.join(args)

    #p = subprocess.Popen( shlex.split(' '.join(args)), stdin=PIPE, stdout=PIPE)
    #p.wait()
    #output = p.communicate("get file.ext")
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
