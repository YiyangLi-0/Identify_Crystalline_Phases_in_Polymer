#!/usr/bin/env python3
import numpy as np
from mpi4py import MPI
# Custom modules.
from . import md_system
from . import bin_table
from . import p2_parameters


def possible_crys_atoms(comm, mds, bins, inp):
    """ Manage worker processes to identify possible crystalline atoms.
    """
    ''' Assign local set of atoms to be handled by current process. '''
    chunked_atoms = md_system.distribute_atoms(range(mds.n_atom), comm.size)
    loc_atoms = chunked_atoms[comm.rank]

    ''' Initialize an numpy array to store the ids of possible crystalline
        atoms (pca) for the current process. '''
    L   = max([len(_) for _ in chunked_atoms])  # Length of pca.
    pca = np.empty(L, dtype = np.int)
    pca.fill(-1)

    ''' Initialize a list to store the contents of pca collected from all
        worker processes. '''
    if comm.rank == 0:
        master_pca = np.empty(comm.size*L, dtype = np.int)
    else:
        master_pca = None

    ''' Identify crystalline atoms and modify pca in current process. '''
    proc = '{}-rank{}'.format(MPI.Get_processor_name(), comm.rank)
    p2_parameters.p2_parameters(mds, bins, inp, loc_atoms, pca, proc)

    ''' Collect pca from worker processes to master process. '''
    comm.Gather(pca, master_pca, root = 0)
    if comm.rank == 0:
        master_pca = list(set(master_pca))  # Trim additional -1s.
        master_pca.remove(-1)
    return master_pca  # 0-indexed list for rank 0, None for other ranks.


if __name__ == '__main__':
    pass
