#!/usr/bin/env python

from matplotlib import pyplot as plt
import sys
import numpy as np

if __name__ == '__main__':
    print("\n\n")
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(100))
    print("            ___________________________           ".center(100))
    print(" __________|  C o d e  b y  Y.J. Choi  |_________ ".center(100))
    print("|______________ ssrokyz@gmail.com _______________|".center(100))
    print("")
    print("*******   This code will show you a histogram   *******".center(100))
    print("useage ==> ./histogram.py 'file' 'colume of data' 'start line' 'end line' 'line interval'                    ".center(100))
    print("    or ==> ./histogram.py 'file' 'colume of data' 'start line' 'end line' 'line interval' 'bin min' 'bin max'".center(100))
    print("note) interval 1 means every line".center(100))
    print("EXAMPLE) ./histogram.py log_gst.txt 3 1 -1 2          ".center(100))
    print("      or ./histogram.py log_gst.txt 3 1 -1 2 -200 -250".center(100))
    print("")
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(100))
    print("")
    if len(sys.argv) is 6 or len(sys.argv) is 8 :
        print(("The Number of arguments(= %d) is correct." %(len(sys.argv)-1)).center(100))
        print("\n")
    else:
        print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(100))
        print("\n")
        sys.exit(1)
    file_name = sys.argv[1]
    column = int(sys.argv[2])
    start_line = int(sys.argv[3])
    end_line = int(sys.argv[4])
    interval = int(sys.argv[5])
    if len(sys.argv) == 8:
        bin_range = True
        bin_min = float(sys.argv[6])
        bin_max = float(sys.argv[7])
    else:
        bin_range = False
        bin_min = None
        bin_max = None
    from ss_util import column2np
    data = column2np(file_name, column, start_line, end_line, interval)
    average = np.average(data)
    std = np.std(data)

    if bin_range:
        n, bins, patches = plt.hist(data, bins=200, range=[bin_min, bin_max], facecolor='purple', alpha=0.70)
    else:
        n, bins, patches = plt.hist(data, bins=200, facecolor='purple', alpha=0.70)
    
    #### plot
    plt.title('histogram')
    plt.xlabel('%d data, average = %.3f, sigma = %.3f' % (len(data), average, std))
    plt.ylabel('population')
    plt.grid(True)
    plt.show()

