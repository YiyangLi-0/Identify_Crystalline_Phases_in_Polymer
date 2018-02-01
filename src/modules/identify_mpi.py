#!/usr/bin/env python
import numpy as np
from mpi4py import MPI
# Custom modules.
from . import md_system
from . import bin_table
from . import p2_parameters


def possible_crys_atoms(comm, mds, bins, inp):
    """ Manage worker processes to identify possible crystalline atoms.
    """
    # Approximate number of local atoms to be handled by a worker process.
    loc_n = mds.n_atom / comm.size

    ''' Generate a list of atom ids to be handled by current worker process. '''
    atoms = range(mds.n_atom)
    if comm.rank < comm.size - 1:
        loc_atoms = atoms[comm.rank*loc_n : (comm.rank+1)*loc_n]
    else:
        loc_atoms = atoms[comm.rank*loc_n :]

    ''' Initialize an numpy array to store the ids of possible crystalline
        atoms (pca) for the current worker process. '''
    L   = max(loc_n, mds.n_atom-(comm.size-1)*loc_n)  # Length of pca.
    pca = np.empty(L, dtype = np.int)
    pca.fill(-1)

    ''' Initialize a list to store the contents of pca collected from all
        worker processes. '''
    if comm.rank == 0:
        master_pca = np.empty(comm.size*L, dtype = np.int)
    else:
        master_pca = None

    ''' Identify crystalline atoms and modify pca in current worker process. '''
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
