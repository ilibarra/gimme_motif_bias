'''
Created on 5/3/2018

DESCRIPTION

@author: ignacio
'''
import threading
from threading import Thread

from multiprocessing import Process
from multiprocessing import Manager

class ThreadingUtils:

    @staticmethod
    def chunks(l, n):
        """Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield l[i:i + n]

    @staticmethod
    def run(function, input_list, n_cores):
        print(('run function with %i cores' % n_cores))
        print(function)
        # print 'with input list of len'
        # print len(input_list)
        # print 'in groups of %d threads' % n_threads

        assert n_cores <= 20
        n_groups = int(len(input_list) / n_cores + 1)
        # print 'n groups', n_groups
        for group_i in range(n_groups):
            start, end = group_i * n_cores, (group_i + 1) * n_cores
            # print 'start', start, 'end', end

            threads = [None] * (end - start)
            for i, pi in enumerate(range(start, min(end, len(input_list)))):
                next_args = input_list[pi]
                # print next_kmer
                threads[i] = Process(target=function, args=next_args)
                # print 'starting process #', i
                threads[i].start()

            # print  threads
            # print 'joining threads...'
            # do some other stuff
            for i in range(len(threads)):
                if threads[i] is None:
                    continue
                threads[i].join()
