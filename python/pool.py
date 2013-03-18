from multiprocessing import Pool, cpu_count

PROCESSES=cpu_count()
_pool = None

def start_pool():
    global _pool
    if _pool is None:
        _pool = Pool(processes=PROCESSES)
    
def stop_pool():
    global _pool
    _pool.close()
    _pool = None

def check_pool():
    return (_pool is not None)

def pool():
    return _pool
