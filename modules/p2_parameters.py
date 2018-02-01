#!/usr/bin/env python
import sys
import numpy as np
# Custom modules.
from . import md_system
from . import bin_table


def p2_parameters(mds, bins, inp, loc_atoms, pca, proc):
    """ Calculate local p2 order parameter for each atom i in loc_atoms (a list
        of local atoms to be handled by a parallel process), and modify pca (a
        local list of possible crystalline atoms).
    """
    for m, i in enumerate(loc_atoms):
        xi = mds.atom_coords[i][:]

        # Loop over all bonded neighbors k of atom i.
        p2i = 0.0  # Local p2 order parameter of atom i.
        n_bond_pair = 0  # Number of the pairs (bondi, bondj)
        for k in mds.bonds[i]:
            xi_neig = mds.atom_coords[k][:]
            bondi = md_system.bond_vector(xi, xi_neig, mds.box_length)

            # Loop over atom j in atom i's surrounding bins.
            for j in bin_table.atom_neighbors(i, mds, bins):
                if j == i:
                    continue
                xj = mds.atom_coords[j][:]
                dij_sqr = md_system.image_distance_sqr(xj, xi, mds.box_length)
                if dij_sqr > inp['cutoff_sqr']:
                    continue
                
                # We consider only one bond containing atom j, to avoid
                # redundant computing. Switching from min() to max() may
                # slightly change the calculated p2i (up to +/- 10%).
                xj_neig = mds.atom_coords[min(mds.bonds[j])][:]
                bondj = md_system.bond_vector(xj, xj_neig, mds.box_length)

                # Aggregating local p2 order parameter p2i.
                cos  = np.dot(bondi, bondj)               
                p2i += 0.5 * (3*cos*cos - 1.0)
                n_bond_pair += 1
    
        if n_bond_pair > 0:
            # Determine if atom i belongs to crystalline phase.
            p2i /= float(n_bond_pair)
            if p2i > inp['p2_thresold']:
                # Modify pca. Valid items in pca will be >= 0.
                pca[m] = i

            # Formated screen output.
            if i % inp['log_interval'] == 0 or i == mds.n_atom - 1:
                is_crys = 'true' if p2i > inp['p2_thresold'] else 'false'
                print '{:>7d} | {:>7.4f}  {:>7}  {:>12}'.format(
                       i, p2i, is_crys, proc)
                sys.stdout.flush()
        else:
            print '  !!! Atom {} has no bond within cutoff {} A'.format(
                     i, inp['cutoff'])
            sys.stdout.flush()


if __name__ == '__main__':
    pass
    
