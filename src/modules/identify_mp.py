#!/usr/bin/env python
import multiprocessing
# Custom modules.
from . import md_system
from . import bin_table
from . import p2_parameters


def possible_crys_atoms(mds, bins, inp):
    """ Manage child processes to identify possible crystalline atoms.
    """
    # Number of child processes.
    n_proc = multiprocessing.cpu_count()
    
    # Approximate number of local atoms to be handled by a child process.
    loc_n = mds.n_atom / n_proc

    ''' Create and submit jobs. '''
    jobs  = [None] * n_proc
    pipes = [None] * n_proc
    atoms = range(mds.n_atom)
    for i in range(n_proc):
        # Generate list of atom ids to be handled by child process i.
        if i < n_proc - 1:
            loc_atoms = atoms[i*loc_n : (i+1)*loc_n]
        else:
            loc_atoms = atoms[i*loc_n :]

        # Initialize a list to store the ids of possible crystalline atoms
        # (pca) for every child process.
        L   = max(loc_n, mds.n_atom-(n_proc-1)*loc_n)  # Length of pca.
        pca = [-1] * L

        # Create and submit jobs.
        parent, child = multiprocessing.Pipe(False)
        p = multiprocessing.Process(
                name   = 'local-rank_{}'.format(i),
                target =  parallel_job,
                args   = (mds, bins, inp, loc_atoms, pca, child))
        p.start()
        jobs[i]  = p
        pipes[i] = parent

    ''' Wait for all child processes to complete. '''
    for p in jobs:
        p.join()

    ''' Collect pca from child processes to master process. '''
    master_pca = []
    for parent in pipes:
        master_pca += parent.recv()     # Receive from the child end of pipe.
    master_pca = list(set(master_pca))  # Trim additional -1s.
    master_pca.remove(-1)
    return master_pca  # 0-indexed.

def parallel_job(mds, bins, inp, loc_atoms, pca, child):
    """ Identify possible crystalline atoms and modify pca.
    """
    proc = multiprocessing.current_process().name
    p2_parameters.p2_parameters(mds, bins, inp, loc_atoms, pca, proc)
    child.send(pca)  # Send to the parent end of pipe.


if __name__ == '__main__':
    pass
