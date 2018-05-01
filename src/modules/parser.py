#!/usr/bin/env python3
"""
  Script for parsing input file.
"""
import os
import unicodedata


def is_number(s):
    """ Test if a string s is a number.
    """
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False

def parser(input_file):
    """ Read from input file to generate a dictionary.
    """
    inp = {'output_dir': os.getcwd() + '/output/'}

    with open(input_file, 'r') as fid:
        while True:
            line = fid.readline()
            if not line:
              break

            line = line.split()
            if len(line) < 2 or '#' in line[0]:
                continue 
            
            key, val = line[:2]
            if is_number(val):
                val = float(val)
                if (val).is_integer():
                    inp[key] = int(val)  # int
                else:
                    inp[key] = val       # float
            else:
                inp[key] = val           # string 

    inp['lmp_data'] = inp['data_dir'] + '/' + inp['lmp_data']
    inp['lmp_trj'] = inp['data_dir'] + '/' + inp['lmp_trj']
    inp['cutoff_sqr'] = inp['cutoff']**2  # angstrom^2
    inp['new_trj_file_w'] = inp['output_dir'] + 'identified_wrapped.lammpstrj'
    inp['new_trj_file_uw'] = inp['output_dir'] + 'identified_unwrapped.lammpstrj'
    inp['crys_atoms_file'] = inp['output_dir'] + 'identified_crystalline_atoms'
    inp['crystallinty_file'] = inp['output_dir'] + 'crystallinty'
    return inp


if __name__ == '__main__':
    pass

