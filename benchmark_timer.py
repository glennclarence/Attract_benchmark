import time
from threading import Lock

class timer:
    def __init__(self, start= 0 ):
        self.start = start
        self.stop = 0
        self.elapsed = 0
    def set_start(self, time_start):
        self.start= time_start
    def set_stop(self, time_stop):
        self.stop = time_stop
        self.elapsed = self.stop - self.start
    def set_stopAppend(self, time_stop):
        self.stop = time_stop
        self.elapsed += self.stop - self.start
    def get_elapsedTime(self):
        return self.elapsed


class measure_benchmark:
    def __init__(self):
        self.timers={}
        self.lock = Lock()
    def timer_add(self, name, start= False):
        self.lock.acquire()
        self.timers[name]=timer()
        if start:
            self.timers[name].set_start( time.time() )
        self.lock.release()

    def timer_addTimer(self, name, timer):
        self.lock.acquire()
        self.timers[name]=timer
        self.lock.release()

    def timer_start(self, name):
        start = time.time()
        self.timers[name].set_start( start )

    def timer_stop(self, name):
        stop = time.time()
        self.lock.acquire()
        self.timers[name].set_stop( stop )
        self.lock.release()
    def timer_stop(self, name):
        self.timers[name].set_stop( time.time() )

    def timer_appendStop(self, name ):
        stop = time.time()
        self.timers[name].set_stopAppend(stop)

    def get_elapsedTime(self, name_timer ):
        return self.timers[name_timer].get_elapsedTime()

    def get_timers(self):
        return self.timers


    def sumup_and_getTimer(self):
        totaltime = 0
        for key, temptimer in self.timers.iteritems():
            totaltime += temptimer.get_elapsedTime()
        newTimer = timer()
        newTimer.elapsed = totaltime
        return newTimer

    def save_benchmark(self, filename_benchmark):
        file = open( filename_benchmark, "w")
        file.write("{:25} \t{}\n\n".format("Timername","Elapsed Time (s)"))
        for name_timer in self.timers:
            line="{:25} \t{:04f}\n".format( name_timer, self.timers[name_timer].get_elapsedTime())
            file.write( line )
        file.close()