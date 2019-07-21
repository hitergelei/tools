#!/usr/bin/env python
import numpy as np

def invert_chem_sym_num(inp):
    """
    inp_list (list or tuple or string or integer) - If type is list(or tuple), the list(or tuple) must contain same typed elements only (chem or num).
    """

    ## Change input type to list
    if type(inp) is list:
        inp_list = inp
    else:
        inp_list = [inp]

    ## Determine dtype
    try:
        inp_arr = np.array(inp_list, dtype=np.int)
    except ValueError:
        print('Chemical symbols are provided'.center(120))
        dtype = 'chem'
        inp_arr = np.array(inp_list, dtype=np.str)
    else:
        print('Atomic numbers are provided'.center(120))
        dtype = 'num'

    ## Invert
    if dtype == 'chem':
        from ase.data import atomic_numbers
        out_list = []
        for i in range(len(inp_arr)):
            out_list.append(atomic_numbers[inp_arr[i]])
        out_arr = np.array(out_list)
    else:
        from ase.data import chemical_symbols
        out_arr = np.array(chemical_symbols)[inp_arr]
    return out_arr

if __name__ == '__main__':
    import sys
    import datetime

    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ POSTECH <<<<<'.center(120))
    print(('code started time: '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('This code will give you the atomic numbers from the chemical symbols'.center(120))
    print('(or vice versa)                                        '.center(120))
    print('')
    print('Useage  ==> ./chemical_symbol_number_inverter.py >chemical symbol or atomic number<'.center(120))
    print('Example 1 ==> ./chemical_symbol_number_inverter.py Ge Sb Te                        '.center(120))
    print('Example 2 ==> ./chemical_symbol_number_inverter.py 32 51 52                        '.center(120))
    print('')
    print('=================================================================================================='.center(120))

    inp_list = []
    for i in range(len(sys.argv)-1):
        inp_list.append(sys.argv[i+1])
    out_arr = invert_chem_sym_num(inp_list)
    print('')
    for i in range(len(inp_list)):
        print('{:2} --> {:2}'.format(inp_list[i], out_arr[i]).center(120))
    print('')
