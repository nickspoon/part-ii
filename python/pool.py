from multiprocessing import Pool, cpu_count
import signal, sys

PROCESSES=cpu_count()
RETIRE_AFTER=None
_pool = None

class InPool:
    def __enter__(self):
        start_pool()
        return pool()
    def __exit__(self, type, value, traceback):
        stop_pool()

def init_worker():
    signal.signal(signal.SIGINT, stop_handler)

def stop_handler(signum, frame):
    raise Exception

def start_pool():
    global _pool
    if _pool is None:
        _pool = Pool(PROCESSES, init_worker, maxtasksperchild=RETIRE_AFTER)
    
def stop_pool():
    global _pool
    if _pool is not None:
        _pool.close()
        _pool = None

def check_pool():
    global _pool
    return (_pool is not None)

def pool():
    global _pool
    return _pool
