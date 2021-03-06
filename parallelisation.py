__author__ = 'Aitor Blanco (aitor.blancomiguez@unitn.it)'
__version__ = '0.01'
__date__ = '12 Aug 2020'

from utils import error
from multiprocessing import Event, Pool

CHUNKSIZE = 1

"""
Places terminating in the global namespace of the worker subprocesses.
This allows the worker function to access `terminating` even though it is
not passed as an argument to the function.
"""
def init_terminating(terminating_):
    global terminating
    terminating = terminating_


"""
Executes each parallelised call of a function

:param arguments: the name of the function and the arguments
:returns: the result of the parallelised execution
"""
def parallel_execution(arguments):
    function, *args = arguments
    if not terminating.is_set():
        try:
            return function(*args)
        except Exception as e:
            terminating.set()
            error("Parallel execution fails: "+str(e), init_new_line=True, exit=True)
    else:
        terminating.set()


"""
Creates a pool for a parallelised function and returns the results of each execution
as a list

:param args: the name of the function and the arguments
:param nprocs: the number of threads to use in the pool
:returns: the result of the parallelised funtion
"""
def execute_pool(args, nprocs):
    terminating = Event()
    with Pool(initializer=init_terminating, initargs=(terminating,), processes=nprocs) as pool:
        try:
            return [_ for _ in pool.imap_unordered(parallel_execution, args, chunksize=CHUNKSIZE)]
        except Exception as e:
            error('Parallel execution fails: '+str(e), init_new_line=True, exit=True)
