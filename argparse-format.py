#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    DESCRIPTION
    """)
    # Positional arguments
    parser.add_argument('inp_file_list', type=str, nargs='+', help='ASE readable atoms list file name. When multiple input files are provided, DOS will averaged.')
    # Optional arguments
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='Image range following python convention. default=":" (e.g.) -n :1000:10')
    return parser.parse_args()

if __name__ == '__main__':
    ## Intro
    import datetime
    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ POSTECH <<<<<'.center(120))
    print(('Code runtime : '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('DESCRIPTION'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()
