#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Read yaml file. Use like ==> "python -i read-yaml.py 'some yaml file'"
    """)
    # Positional arguments
    parser.add_argument('yaml_file', type=str, help='yaml file.')

    return parser.parse_args()

if __name__ == '__main__':
    ## Intro
    import datetime
    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ Phys. Dep. of POSTECH in Korea <<<<<'.center(120))
    print(('Code runtime : '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('Read yaml file. Use like ==> "python -i read-yaml.py some_yaml_file"'.center(120))
    print('                data: yaml data     ')
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args = argparse()
    yaml_file = args.yaml_file

    import yaml
    try:
        from yaml import CLoader as Loader
    except ImportError:
        from yaml import Loader
    with open(yaml_file) as f:
        data = yaml.load(f.read(), Loader=Loader)
