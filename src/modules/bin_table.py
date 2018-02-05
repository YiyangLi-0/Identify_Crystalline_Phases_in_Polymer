#!/usr/bin/env python


class bintable:
    """ Class to generate a bin table for a MD system mds.
    """
    def __init__(self):
        self.xbins = 0   # Number of bins in x direction
        self.ybins = 0   # Number of bins in y direction
        self.zbins = 0   # Number of bins in z direction
        self.bins  = []  # Nested list, each sublist contains all the atom-ids
                         # belonging to a bin.

    def find_bin(self, mds, i):
        """ Find the bin-id that atom i belongs to.
        """
        r  = mds.atom_coords[i]  # Wrapped coordinates.
        
        # Index of bin in x direction.
        bx = (r[0] - min(mds.box[0])) // (mds.box_length[0]/float(self.xbins))
        bx = min(int(bx), self.xbins-1)

        # Index of bin in y direction.
        by = (r[1] - min(mds.box[1])) // (mds.box_length[1]/float(self.ybins))
        by = min(int(by), self.ybins-1)

        # Index of bin in z direction.
        bz = (r[2] - min(mds.box[2])) // (mds.box_length[2]/float(self.zbins))
        bz = min(int(bz), self.zbins-1)
        return bx + by * self.xbins + bz * self.xbins * self.ybins


def initialize_bins(mds, cutoff):
    """ Construct a bin table for a MD system mds, using the specified
        cutoff distance.
    """
    bt = bintable()

    ''' Determine the number of bins in each direction. '''
    bt.xbins = max( int(mds.box_length[0] // cutoff), 1 )
    bt.ybins = max( int(mds.box_length[1] // cutoff), 1 )
    bt.zbins = max( int(mds.box_length[2] // cutoff), 1 )

    ''' Determine which atom belongs to which bin. '''
    nbin = bt.xbins * bt.ybins * bt.zbins
    bt.bins = [[] for _ in range(nbin)]
    for i in range(mds.n_atom):  # Atom-ids are 0-indexed.
        bin_id = bt.find_bin(mds, i)
        bt.bins[bin_id].append(i)
    return bt

def atom_neighbors(m, mds, bt):
    """ Return a list of atoms that are within the nearest surrounding bins
        to atom m.
    """
    ''' Get the indices (bx, by, bz) of the bin that contains atom m. '''
    bin_id = bt.find_bin(mds, m)
    bx =  bin_id % bt.xbins
    by = (bin_id // bt.xbins) % bt.ybins
    bz =  bin_id // bt.xbins // bt.ybins

    ''' Get the indices (nx, ny, nz) of nearest bins surrounding the bin
        (bx, by, bz), then take the atom-ids in these bins. '''
    atoms = []
    for i in [-1, 0, 1]:
        if i == 1 and bt.xbins < 3:
            continue
        nx = (bx + i + bt.xbins) % bt.xbins
        for j in [-1, 0, 1]:
            if j == 1 and bt.ybins < 3:
                continue
            ny = (by +j + bt.ybins) % bt.ybins
            for k in [-1, 0, 1]:
                if k == 1 and bt.zbins < 3:
                    continue
                nz = (bz + k + bt.zbins) % bt.zbins
                bin_id = nx + ny * bt.xbins + nz * bt.xbins * bt.ybins
                atoms.append(bt.bins[bin_id])

    # Flatten 2D list.
    atoms = [a for atoms_sublist in atoms for a in atoms_sublist]
    return atoms


if __name__ == '__main__':
    pass
    