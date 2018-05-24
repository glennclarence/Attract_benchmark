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

    #def run(self):
        #while True:
            #if not  self.queue.full():
                #item = random.randint(1, 10)
                #self.queue.put(item)
                #logging.debug('Putting ' + str(item)
                 #             + ' : ' + str(q.qsize()) + ' items in queue')
                #time.sleep(0.08)
        return

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
        return

    def run(self):
        while ( self.signal or not  self.queue.empty() ) and not self.kill:
            if not  self.queue.empty():
                item =  self.queue.get()
                self.target( item, self.name, self.args )
               # logging.debug(self.name + ' is Getting ' + str(item)
               #              + ' : ' + str( self.queue.qsize()) + ' items in queue')
               # time.sleep(random.random())
        return

    def kill(self):
        self.kill = True

    def stop_if_done(self):
        self.Signal = False

def compute( queue_item, threadId, args = None):
        string = "thread {} item {}".format(threadId, queue_item)
        print string
        #logging.debug(self.name + ' is Getting ' + str(item)
        #              + ' : ' + str(queue.qsize()) + ' items in queue')
        #time.sleep(random.random())



#
# if __name__ == '__main__':
#     p = ProducerThread(quq, name='producer')
#     c = []
#
#     for i in range(8):
#         c.append( ConsumerThread( quq, name='consumer'+str(i) , target = compute, args=quq) )
#
#     p.start()
#     time.sleep(2)
#     for i in range(5):
#         c[i].start()
#     time.sleep(2)
#     p.add( 1 )
#     p.add(2)
#     p.add(3)
#     p.add(4 )
#     p.add(5)
#     p.add(6)
#     p.add(7)
#     p.add(8 )


