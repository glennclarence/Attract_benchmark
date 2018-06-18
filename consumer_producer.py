import threading
import time
import logging
import random
import Queue

logging.basicConfig(level=logging.DEBUG,
                    format='(%(threadName)-9s) %(message)s', )

BUF_SIZE = 10
quq = Queue.Queue(BUF_SIZE)


class ProducerThread(threading.Thread):
    def __init__(self, queue, group=None, target=None, name=None,
                 args=(), kwargs=None, verbose=None):
        super(ProducerThread, self).__init__()
        self.target = target
        self.name = name
        self.queue = queue
        self.signal = True


    # def run(self):
    #     while self.Signal is True:
    #         if not  self.queue.full():
    #             if not self.feeder.empty():
    #                 item = self.feeder.get()
    #                 self.queue.put(item)
    #             #logging.debug('Putting ' + str(item)
    #              #             + ' : ' + str(q.qsize()) + ' items in queue')
    #             #time.sleep(0.08)
    #     return

    def add(self, item):
        self.queue.put(item)


class ConsumerThread(threading.Thread):
    def __init__(self, queue, group=None, target=None, name=None,
                 args=(), kwargs=None, verbose=None):
        super(ConsumerThread, self).__init__()
        self.target = target
        self.args = args
        self.name = name
        self.queue = queue
        self.signal = True
        self.kill = False
        self.computing = False
        return

    def run(self):
        while ( self.signal or not  self.queue.empty() ) and not self.kill:
            if not  self.queue.empty():
                self.computing = True
                item =  self.queue.get()
                self.target( item, self.name, self.args )
                self.computing = False
               # logging.debug(self.name + ' is Getting ' + str(item)
               #              + ' : ' + str( self.queue.qsize()) + ' items in queue')
               # time.sleep(random.random())
        return

    def kill(self):
        self.kill = True

    def stop_if_done(self):
        self.signal = False

    def is_done(self):
        return not (self.computing or self.queue.qsize() > 0 or self.signal  )

