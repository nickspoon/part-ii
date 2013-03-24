from multiprocessing import Pool, cpu_count
import signal, sys

PROCESSES=cpu_count()
_pool = None

def init_worker():
    signal.signal(signal.SIGINT, stop_handler)

def stop_handler(signum, frame):
    raise Exception

def start_pool():
    global _pool
    if _pool is None:
        _pool = Pool(PROCESSES, init_worker)
    
def stop_pool():
    global _pool
    _pool.close()
    _pool = None

def check_pool():
    return (_pool is not None)

def pool():
    if _pool is None:
        start_pool()
    return _pool
