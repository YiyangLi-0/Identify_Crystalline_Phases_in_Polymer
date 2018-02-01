#!/usr/bin/env python
import numpy as np
import copy


class md_system:
    """ Class to store molecular information of a linear polyethylene system.
        Atom-ids are 0-indexed, i.e., self.atom_types[i] is the type of the
        i-th atom, and self.bonds[i] stores all atoms connected to atom i in
        a list. Here 'atom' means 'coarse-grained bead', and simulation box
        is orthorhombic.
    """
    def __init__(self):
        self.n_atom       = 0   # Number of atoms in MD system.
        self.atom_coords  = []  # Wrapped coordinate of atoms.
        self.atom_types   = []  # Atom types of atoms.
        self.is_crys_atom = []  # If atoms are identified as belonging
                                # to crystalline phase, 0-false, 1-true.
        self.n_chain      = 0   # Number of chains in MD system.
        self.chain        = []  # Which chain the atoms belongs to.
        self.chain_ids    = []  # The ids of chains in MD system, read from
                                # LAMMPS data file, may not start from 0.
        self.chains       = []  # List containing the atoms belongs to each
                                # linear chain, atoms are sorted according to
                                # atom connectivity.
        self.n_bond       = 0   # Number of bonds in MD system.
        self.bonds        = []  # Which atoms the current atom connects with.
        self.box          = []  # Bounds of simulation box in each direction.
                                # Only orthorhombic box is considered here.
        self.box_length   = []  # Lenght of box in each direction.
        self.time         = 0   # Time read from LAMMPS trajectory file.

    def crystalline_atoms(self):
        """ Return a list of all crystalline atoms in MD system.
        """
        return [i for i in range(self.n_atom) if self.is_crys_atom[i]]

    def chain_ends(self, c):
        """ Rerurn the end-atoms of a linear chain with chain-id c.
        """
        atoms = [i for i, x in enumerate(self.chain) if x == c]
        ends  = []
        for atom in atoms:
            if len(self.bonds[atom]) == 1:
                ends.append(atom)
        return ends

    def atoms_in_chain(self, c):
        """ Rerurn all the atoms belonging to a chain with chain-id c.
            Atoms are sorted according to their connectivity.
        """
        ends  = self.chain_ends(c)
        atoms = [ends[0]]               # Start from a random end atom
        queue = self.bonds[ends[0]][:]  # Avoid modifying self.bonds
        while queue:
            next = queue[-1]
            queue.pop(-1)
            atoms.append(next)
            for neighbor in self.bonds[next]:
                if neighbor not in atoms:
                    queue.append(neighbor)
        return atoms

    def find_all_segments(self, atoms):
        """ Identify all atom-segments (due to atom connectivity) from the
            input list of atoms. No duplicate atom-ids should appear in atoms,
            and note that atoms will be modified.
        """
        segments = []
        while atoms:
            segment = self.find_one_segment(atoms)
            segments.append(segment)
            # Remove segemnt from atoms to avoid duplicate work.
            atoms = list(set(atoms)-set(segment))
        return segments
    
    def find_one_segment(self, atoms):
        """ Identify one atom-segment from the input list of atoms.
        """
        ini     = atoms[0]
        segment = [ini]
        queue   = self.bonds[ini][:]
        while queue:
            next = queue[-1]
            queue.pop(-1)
            if (next not in atoms) or (next in segment):
                continue
            segment.append(next)
            for neighbor in self.bonds[next]:
                if neighbor in atoms:
                    if neighbor not in segment:
                        queue.append(neighbor)
        return segment


def read_md_system(lmp_data):
    """ Read MD system information from LAMMPS data file.
    """
    mds = md_system()
    with open(lmp_data, "r") as fid:
        ''' Read number of atoms and number of bonds in MD system. '''
        while True:
            line = fid.readline()
            if 'xlo' in line:
                break
            if 'atoms' in line:
                mds.n_atom = int(line.split()[0])
            if 'bonds' in line:
                mds.n_bond = int(line.split()[0])

        ''' Initialize mds.is_crys_atom for later use. '''
        mds.is_crys_atom = [0] * mds.n_atom

        ''' Read box bounds. '''
        for i in range(3):
            if i == 0:
                # x direction
                line = line.split()
            else:
                # y, z directions
                line = fid.readline().split()
            mds.box.append( [float(j) for j in line[0:2]] )
            mds.box_length.append( abs(mds.box[i][1]-mds.box[i][0]) )
    
        ''' Read atom types and coordinates. Each line follows this format:
            atom-id, mol-id, atom-type, x, y, z '''
        # Fast forward.
        while True:
            line = fid.readline()
            if line.startswith('Atoms'):
                fid.readline()
                break
        # Read types and coordinates.
        mds.atom_coords = np.zeros((mds.n_atom, 3))
        mds.atom_types  = [-1] * mds.n_atom
        mds.chain       = [-1] * mds.n_atom
        for i in range(mds.n_atom):
            line = fid.readline().split()
            atom_id = int(line[0]) - 1  # Convert ID to 0-indexed.
            mds.chain[atom_id] = int(line[1])
            mds.atom_types[atom_id] = int(line[2])
            x = [float(s) for s in line[3:6]]
            wrap_atom_into_box(x, mds.box, mds.box_length)
            mds.atom_coords[atom_id] = x

        ''' Read bond information. Each line follows this format:
            bond-id, bond-type, atom-id1, atom-id2 '''
        # Fast forward.
        while True:
            line = fid.readline()
            if line.startswith('Bonds'):
                fid.readline()
                break
        # Read and generate bond table.
        mds.bonds = [set() for _ in range(mds.n_atom)]
        for _ in range(mds.n_bond):
            line = fid.readline().split()
            i, j = [int(s)-1 for s in line[2:4]]  # Convert IDs to 0-indexed.
            mds.bonds[i].add(j)
            mds.bonds[j].add(i)
        mds.bonds = [sorted(list(b)) for b in mds.bonds]

        ''' Call established mds.bonds to identify chains in MD system. '''
        mds.chain_ids = sorted(list(set(mds.chain)))
        for c in mds.chain_ids:
            mds.chains.append(mds.atoms_in_chain(c))
        mds.n_chain = len(mds.chains)
        # Check consistency.
        assert  mds.n_chain == len(mds.find_all_segments(range(mds.n_atom))), \
                "Atom connectivity is inconsistent in {}".format(lmp_data)
    return mds

def wrap_atom_into_box(x, box, box_length):
    """ Put atom x inside orthorhombic box using periodic boundary conditions.
    """
    for i in range(3):
        while x[i] < min(box[i]):
            x[i] += box_length[i]
        while x[i] > max(box[i]):
            x[i] -= box_length[i]

def nearest_atom_image(x, x_ref, box_length):
    """ Return the image coordinate of atom x in the nearest vicinity of
        reference atom x_ref, according to PBC. No modification is made
        to x or x_ref.
    """
    dl = np.array(x) - np.array(x_ref)
    for i in range(3):
        while dl[i] < -0.5 * box_length[i]:
            dl[i] += box_length[i]
        while dl[i] > 0.5 * box_length[i]:
            dl[i] -= box_length[i]
    return x_ref + dl

def image_distance_sqr(x1, x2, box_length):
    """ Return the squared image distance between atoms x1 and x2.
    """
    x2_img = nearest_atom_image(x2, x1, box_length)
    bond = x2_img - x1
    return np.dot(bond, bond)

def bond_vector(xi, xi_neig, box_length):
    """ Determine the normalized bond vector between atoms i and j.
    """
    bond  = xi - nearest_atom_image(xi_neig, xi, box_length)
    bond /= np.linalg.norm(bond)
    return bond

def read_frame_positions(lmp_trj):
    """ Read stream positions in trajectory file corresponding to
        time-step and atom-data.
    """
    ts_pos, data_pos = [], []
    with open(lmp_trj, 'r') as fid:
        while True:
            line = fid.readline()
            if not line:
                break
            if line.startswith('ITEM: TIMESTEP'):
                ts_pos.append(fid.tell())
            elif line.startswith('ITEM: ATOMS id'):
                data_pos.append(fid.tell())
    return ts_pos, data_pos

def update_md_system(mds_ref, lmp_trj, t_pos, d_pos, dt):
    """ Return a new instance of md_system class by reading from LAMMPS
        trajectory file at the specified location (t_pos, d_pos).
    """
    mds = copy.deepcopy(mds_ref)
    with open(lmp_trj, 'r') as fid:
        ''' Identify the column indices of 'x', 'y', 'z' in trajectory file. '''
        scaled_coordinates = False
        x_str, y_str, z_str = 'x', 'y', 'z'
        while True:
            line = fid.readline()
            if line.startswith('ITEM: ATOMS id'):
                if 'xs' in line:
                    scaled_coordinates = True
                    x_str, y_str, z_str = 'xs', 'ys', 'zs'
                line = line.split()
                xid = line.index(x_str) - 2
                yid = line.index(y_str) - 2
                zid = line.index(z_str) - 2
                break

        ''' Update simulation time and box bounds. '''
        fid.seek(t_pos)
        mds.time = float(fid.readline()) * dt
        for i in range(3):
            fid.readline()
        for i in range(3):
            line = fid.readline().split()
            mds.box[i] = [float(j) for j in line[0:2]]
            mds.box_length[i] = abs(mds.box[i][1] - mds.box[i][0])

        ''' RE-initialize mds.is_crys_atom. '''
        mds.is_crys_atom = [0] * mds.n_atom

        ''' Update atom coordinates. '''
        fid.seek(d_pos)
        while True:        
            line = fid.readline()
            if not line or line.startswith('ITEM: TIMESTEP'):
                break
            atom_id = int(line.split()[0]) - 1  # Convert to 0-indexed
            xs = [float(p) for p in line.split()[xid : zid+1]]
            if scaled_coordinates:
                # If the coordinates xs in trajectory file are scaled between
                # [0., 1.], convert them to actual values.
                x = mds.box[0][0] + xs[0] * mds.box_length[0]
                y = mds.box[1][0] + xs[1] * mds.box_length[1]
                z = mds.box[2][0] + xs[2] * mds.box_length[2]
                xs = [x,y,z]
            wrap_atom_into_box(xs, mds.box, mds.box_length)
            mds.atom_coords[atom_id] = xs
    return mds

def write_trj(mds, trj_target):
    """ Write a LAMMPS trajectory showing which phase (amorphous or
        crystalline) every atom belongs to.
    """
    with open(trj_target, 'a') as wid:
        wid.write('ITEM: TIMESTEP\n')
        wid.write('{}\n'.format(mds.time))

        wid.write('ITEM: NUMBER OF ATOMS\n')
        wid.write('{}\n'.format(mds.n_atom))

        wid.write('ITEM: BOX BOUNDS pp pp pp\n')
        for i in range(3):
            wid.write('{} {}\n'.format(mds.box[i][0], mds.box[i][1]))

        ''' VMD reads atom_type only at the first time step in trajectory file,
            so atom_types keep constants for later time steps. To make sure
            different atoms can be rendered with corresponding colors
            in VMD, I use the user field 'vx' to store dynamically changing
            atom_types. '''
        wid.write('ITEM: ATOMS id mol type x y z vx\n')
        for i in range(mds.n_atom):
            atom_id    = i + 1            # convert to 1-indexed.
            molecule   = mds.chain[i] + 1 # convert to 1-indexed.
            atom_type  = mds.atom_types[i]
            x,y,z      = mds.atom_coords[i]
            is_crystal = mds.is_crys_atom[i]
            row = (atom_id, molecule, atom_type, x, y, z, is_crystal)
            wid.write(' %6d %4d %2d %.6f %.6f %.6f %d \n' %row)


if __name__ == '__main__':
    pass
