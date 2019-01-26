import numpy as np

if __name__ == '__main__':
    import sys
    import datetime

    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print("")
    print(">>>>> Code by Young Jae Choi @ POSTECH <<<<<".center(120))
    print(("code started time: "+time).center(120))
    print("")
    print("==================================================================================================".center(120))
    print("")
    print("Useage  ==> python -i pckl-read.py 'pckl file'".center(120))
    print("Example ==> python -i pckl-read.py phonon.pckl".center(120))
    print('')
    print('>>>> Pre-definition <<<<'.center(120))
    print('     pckl    : Pickle class              '.center(120))
    print('     obj     : Loaded pickle class object'.center(120))
    print('     keys    : obj.__dict__.keys()       '.center(120))
    print("")
    print("==================================================================================================".center(120))
    
    if len(sys.argv) == 2:
        pckl_file = sys.argv[1]
    else:
        raise ValueError("*****ERROR***** The number of arguments is not correct *****ERROR*****")

    if sys.version_info[0] == 3:
        import pickle as pckl
    else:
        import cPickle as pckl
    with open(pckl_file, 'rb') as f:
        obj = pckl.load(f)
    keys = obj.__dict__.keys()
