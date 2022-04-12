__author__ = ('Leonard Dubois (leonard.dubois@unitn.it), '
              'Aitor Blanco (aitor.blancomiguez@unitn.it)')
__version__ = '0.01'
__date__ = '27 Aug 2020'

"""
Downloading UniRef diamond databases for Uniref Pangenome builder
"""

import os, subprocess, sys, time, bz2
import argparse as ap
from urllib.request import urlretrieve, urlcleanup
from utils import info


NR90_DMND = ["UniRef90_201906.dmnd", "http://cmprod1.cibio.unitn.it/databases/PanPhlAn/UniRef90_201906.dmnd"]
NR50_DMND = ["UniRef50_201906.dmnd", "http://cmprod1.cibio.unitn.it/databases/PanPhlAn/UniRef50_201906.dmnd"]
MAP_90_TO_50 = ["UniRef90to50_201906.tsv", "http://cmprod1.cibio.unitn.it/databases/PanPhlAn/UniRef90to50_201906.tsv"]

# ------------------------------------------------------------------------------
"""
Reads and parses the command line arguments of the script.
:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('-o', '--output', type = str, default = ".",
                    help='')
    p.add_argument('-v', '--verbose', action='store_true',
                    help='Show progress information')
    return p.parse_args()


def check_args(args):
    if not os.path.exists(args.output):
        os.mkdir(args.output)

# ------------------------------------------------------------------------------
#   DOWNLOAD FILES
# ------------------------------------------------------------------------------
def byte_to_megabyte(byte):
    """Convert byte value to megabyte
    """
    return (byte / 1048576)
    
class ReportHook():

    def __init__(self):
        self.start_time = time.time()

    def report(self, blocknum, block_size, total_size):
        """Print download progress message
        """
        if blocknum == 0:
            self.start_time = time.time()

            if total_size > 0:
                info("Downloading file of size: {:.2f} MB\n".format(byte_to_megabyte(total_size)))
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

            if total_size > 0:
                percent_downloaded = total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stderr to overwrite stderr
                download_rate = total_downloaded / (time.time() - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += ("{:3.2f} %  {:5.2f} MB/sec {:2.0f} min {:2.0f} sec "
                           .format(percent_downloaded, byte_to_megabyte(download_rate),
                                   estimated_minutes, estimated_seconds))

            status += "        \r"
            info(status)


def download(url, download_file, overwrite=False, verbose=False):
    """Download a file from a url
    """
    urlcleanup()
    if (not os.path.isfile(download_file)) or overwrite:
        try:
            if verbose:
                info('Downloading "{}" to "{}"\n'.format(url, download_file))
            urlretrieve(url, download_file, reporthook=ReportHook().report)
            info('\n')
        except EnvironmentError:
            info('unable to download "{}"'.format(url), exit=True)
    else:
        info('File "{}" already present\n'.format(download_file))



# ------------------------------------------------------------------------------
#   MAIN
# ------------------------------------------------------------------------------
def main():
    if not sys.version_info.major == 3:
        sys.stderr.write('[E] Python version: ' + sys.version)
        sys.exit('[E] This software uses Python 3, please update Python')
    args = read_params()
    check_args(args)
    
    for filename, url in [NR90_DMND, NR50_DMND, MAP_90_TO_50] :
        download(url, os.path.join(args.output, filename), verbose =  args.verbose )
    

if __name__ == '__main__':
    start_time = time.time()
    main()
    mins_elapsed = round((time.time() - start_time) / 60.0, 2)
    print('[TERMINATING...] ' + __file__ + ', ' + str(mins_elapsed) + ' minutes.')
