#!/usr/bin/env python
""" 
  Identify crystalline atoms and crystalline chunks in a trajectory of MD
  system. Here 'atom' means 'coarse-grained bead'. Distributed computing is
  separately implemented in MPI mode and multiprocessing mode. To further
  improve efficiency, the bin-table technique is used in both modes.
  
  -----------------------------------------------
  MPI mode:
  
  Run distributed computing on all the worker nodes in your cluster network.
  The master node can also be utilized along side the worker nodes.
  (This can be changed by modifying the 'host' file.) Usage:
  
    cd run
    mpiexec -np <n> python ../src/main.py mpi
  or
    cd run
    mpiexec --hostfile ./hosts -np <n> python ../src/main.py mpi
    
  where <n> is the number of physical cores in your cluster. To use the
  '--hostfile' option, you need to modify the 'host' file to accommodate
  the setting of your cluster.

  -----------------------------------------------
  Multiprocessing mode:
  
  Run parallel computing only on the master node. Usage:
  
    python ../src/main.py mp

"""
import os, sys
# Custom modules
from modules import parser
from modules import md_system
from modules import bin_table
from modules import postprocess


def main():
    if sys.argv[-1] == 'mpi':
        mpi_routine('input')
        
    elif sys.argv[-1] == 'mp':
        mp_routine('input')
        
    else:
        print '!! Wrong parameter provided, should be either \"mpi\" or \"mp\"'
        exit()

def mpi_routine(input):
    """ Using MPI.
        Unless otherwise stated, codes are simultaneously executed on all the
        cores of the nodes listed in the 'hosts' file.
    """
    from mpi4py import MPI
    from modules import identify_mpi  # Custom module

    # Initialize MPI communicator.
    comm = MPI.COMM_WORLD

    ''' Generate input and preprocess output folder. '''
    inp = parser.parser('input')
    if comm.rank == 0:
        check_output_dir(inp)

    ''' Generate reference MD system from LAMMPS data file. '''
    mds_ref = md_system.read_md_system(inp['lmp_data'])
    if comm.rank == 0:
        print 'LAMMPS data file: {}'.format(inp['lmp_data'])

    ''' Identify crystalline atoms at every time step in trajectory file. '''
    ts, ds = md_system.read_frame_positions(inp['lmp_trj'])
    if comm.rank == 0:
        print 'LAMMPS trajectory file: {}\n'.format(inp['lmp_trj'])

    ''' Loop over time steps. '''
    for k in range(len(ts)):
        # Generate MD system and bin table at the current time step.
        comm.Barrier()
        mds = md_system.update_md_system(
                  mds_ref, inp['lmp_trj'], ts[k], ds[k], inp['dt'])
        bins = bin_table.initialize_bins(mds, inp['cutoff'])
        if comm.rank == 0:
            print_header(mds.time)

        # Identify possible crystalline atoms (pca).
        pca = identify_mpi.possible_crys_atoms(comm, mds, bins, inp)

        # Postprocess on the master process.
        if comm.rank == 0:
            postprocess.postprocess(mds, inp, pca)

def mp_routine(input):
    """ Using multiprocessing.
        Unless otherwise stated, codes are executed only on one core of the
        master node.
    """
    from modules import identify_mp  # Custom module

    ''' Generate input and preprocess output folder. '''
    inp = parser.parser(input)
    check_output_dir(inp)

    ''' Generate reference MD system from LAMMPS data file. '''
    mds_ref = md_system.read_md_system(inp['lmp_data'])
    print 'LAMMPS data file: {}'.format(inp['lmp_data'])

    ''' Identify crystalline atoms at every time step in trajectory file. '''
    ts, ds = md_system.read_frame_positions(inp['lmp_trj'])
    print 'LAMMPS trajectory file: {}\n'.format(inp['lmp_trj'])

    ''' Loop over time steps. '''
    for k in range(len(ts)):
        # Generate MD system and bin table at current time step.
        mds = md_system.update_md_system(
                  mds_ref, inp['lmp_trj'], ts[k], ds[k], inp['dt'])
        bins = bin_table.initialize_bins(mds, inp['cutoff'])
        print_header(mds.time)

        # Parallelly identify possible crystalline atoms (pca).
        pca = identify_mp.possible_crys_atoms(mds, bins, inp)

        # Postprocess on the master proces.
        postprocess.postprocess(mds, inp, pca)

def check_output_dir(inp):
    """ Preprocessing output folder and its content.
    """
    if not os.path.exists(inp['output_dir']):
        os.mkdir(inp['output_dir'])
        
    for target in [inp['new_trj_file'], 
                   inp['crys_atoms_file'],
                   inp['crystallinty_file']]:
        if os.path.exists(target):
            os.remove(target)

def print_header(time):
    """ Print header of table.
    """
    print 'Time: {} ns'.format(time)
    print '{:>7} | {:>7}  {:>7}  {:>12}'.format('atom_id','p2','is_crys','proc')
    print '{:->7} | {:->7}  {:->7}  {:->12}'.format('', '', '', '')
    sys.stdout.flush()


if __name__ == '__main__':
    main()
