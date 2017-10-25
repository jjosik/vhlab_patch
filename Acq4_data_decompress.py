#-------------------------------------------------------------------------------
# Name:        Acq4_data_decompress
# Purpose:
#
# Author:      Jason Osik, vhlab
#
# Created:     01/02/2015
# Copyright:   (c) vhlab 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------


def main():
    import os
    import pyqtgraph.metaarray as ma
    for (dirname, dirshere, fileshere) in os.walk(r'/Users/osik_mac/Documents/MATLAB/data_spontReorg/1D/2015.3.27_s0c5/CCIV-Long_004'):
        for filename in fileshere:
            if filename.endswith('.ma'):
                os.chdir(dirname)
                data = ma.MetaArray(file=filename)
                new_filename = "Clamp1_uncomp.ma"
                data.write(new_filename, compression=None)

if __name__ == '__main__':
    main()
