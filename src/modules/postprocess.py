#!/usr/bin/env python
import sys
# Custom modules.
from . import md_system


def postprocess(mds, inp, pca):
    """ Complete the identification of crystalline atoms from possible
        crystalline atoms (pca).
    """
    ''' Filter crystalline_atoms due to chunk length. '''
    sys.stdout.flush()
    print 'Filtering crystalline atoms...'
    chunks = mds.find_all_segments(pca)
    chunks = [chunk for chunk in chunks
              if len(chunk) >= inp['min_chk_length']]

    with open(inp['crys_atoms_file'], 'a') as wid:
        wid.write('\n# time: {} ns\n'.format(mds.time))
        wid.write('Atom ids are 0-indexed\n')
        for chunk in chunks:
            for i in chunk:
                wid.write('{} '.format(i))
            wid.write('\n')

    ''' Generate mds.is_crys_atom and write two new lammps trajectory files.
    '''
    mds.is_crys_atom = [0] * mds.n_atom
    for i in range(mds.n_atom):
        if any(i in chunk for chunk in chunks):
            mds.is_crys_atom[i] = 1

    md_system.write_trj(mds, inp['new_trj_file_w'])   # Wrapped.
    md_system.write_trj(mds, inp['new_trj_file_uw'])  # UnWrapped.
    sys.stdout.flush()
    print 'Done'

    ''' Compute crystallinity of identified MD system. '''
    n_ca = sum([len(c) for c in chunks])
    crystallinity = float(n_ca) / mds.n_atom
    sys.stdout.flush()
    print 'Crystalline atoms:  {} / {}'.format(n_ca, mds.n_atom)
    print 'Crystallinty: {:>4.3}\n'.format(crystallinity)

    with open(inp['crystallinty_file'], 'a') as wid:
        wid.write('{:>6.2f} ns: {:>6.3f}\n'.format(
                   mds.time, crystallinity))


if __name__ == '__main__':
    pass
