# License: GNU Affero General Public License v3 or later

import os
import time
import traceback

from multiprocessing import Process,Queue
from log import return_logger


logger = return_logger(__name__, False)

def put_jobs(jobs,queue,nr):
    for job in jobs[0:nr]:
        queue.put(job)
    return jobs[nr:]
    
def clear_queue(queue):
    while queue.qsize() > 0:
        _ = queue.get()
    queue.close()
    
def stop_processes(processes):
    for process in processes:
        process.terminate()
        process.join()

def operator(jobs,nr_workers,target_function,target_args,quiet=False,sleeptime=1,multiplier=5):
    
    job_queue = Queue()
    results_queue = Queue(maxsize=nr_workers)
    args = [target_args,job_queue,results_queue,target_function,quiet]
    jobs = put_jobs(jobs,job_queue,multiplier*nr_workers)
    
    workers = []
    for i in range(nr_workers):
        worker = Process(target=worker_wrapper,args=args)
        workers.append(worker)
        
    for w in workers:
        w.start()
    results = []
    put_ends = False
    counter_print = 0
    print_status = False
    multiplier_to_print = 100
    while any([w.is_alive() for w in workers]):
        if print_status:
            logger.debug('Working with target function: %s' %target_function.func_name)
            logger.debug('%i workers still alive' %(len([w.is_alive() for w in workers])))
            logger.debug('Work queue: %i; all_groups: %i' %(job_queue.qsize(),len(jobs)))
        while results_queue.qsize() > 0:
            res = results_queue.get()
            if res == None:
                break
            elif type(res) == ValueError:
                logger.error('Encountered error - terminating')
                stop_processes(workers)
                clear_queue(job_queue)
                raise res
            else:
                results.append(res)
                if print_status:
                    logger.debug('New results. Remaining in q: %i' %results_queue.qsize())
        if print_status:
            logger.debug('Continuing')
        if job_queue.qsize() < nr_workers:
            if len(jobs) > 0:
                jobs = put_jobs(jobs,job_queue,multiplier*nr_workers)
            elif not put_ends:
                for _ in range(nr_workers):
                    job_queue.put(False)
                put_ends = True
        # Print only once every 10 times it checks the queue
        counter_print += 1
        if counter_print % multiplier_to_print == 0:
            print_status = True
        else:
            print_status = False
        time.sleep(sleeptime)
    while not results_queue.empty():
        res = results_queue.get()
        results.append(res)
    logger.debug('Joining workers')
    for worker in workers:
        worker.join()
    logger.debug('Final run')
    while not results_queue.empty():
        res = results_queue.get()
        results.append(res)
    return results
    
def worker_wrapper(args,job_queue,results_queue,function,quiet):
    process_id = os.getpid()
    if not quiet:
        logger.debug('Starting process id: %s'  %os.getpid())
    while True:
        try:
            job = job_queue.get()
        except:
            time.sleep(10)
            continue
        if job == False:
            break
        if type(job) == str:
            logger.debug('Process id %s starting work on job %s' %(process_id, job))
        elif type(job) in (tuple,list):
            logger.debug('Process id %s starting work on job starting with %s' %(process_id, job[0]))
        elif type(job) == dict:
            logger.debug('Process id %s starting work on job starting with %s' %(process_id, job.keys()[0]))
            
        job_args = [job] + args
        try:
            result = function(*job_args)
            results_queue.put(result)
        except:
            error = traceback.format_exc()
            logger.error('Process id %s ran into an error' %process_id)
            results_queue.put(ValueError('Process id %s ran into an error:\n%s' %(process_id, error)))
    logger.debug('Process %s exiting' %os.getpid())

