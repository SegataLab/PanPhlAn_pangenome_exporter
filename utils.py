__author__ = 'Aitor Blanco (aitor.blancomiguez@unitn.it)'
__version__ = '0.01'
__date__ = '12 Aug 2020'

import time, sys

def info(s, init_new_line=True, exit=False, exit_value=0, time_stamp=True):
    if init_new_line:
        sys.stdout.write('\n')
    if time_stamp:
        timestamp = time.ctime(int(time.time()))
        sys.stdout.write('{} '.format(timestamp))
    sys.stdout.write('{}'.format(s))
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)


def error(s, init_new_line=False, exit=False, exit_value=1, time_stamp=True):
    if init_new_line:
        sys.stderr.write('\n')
    if time_stamp:
        timestamp = time.ctime(int(time.time()))
        sys.stderr.write('{} '.format(timestamp))
        
    sys.stderr.write('[e] {}\n'.format(s))
    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)


def get_genome_name(genome):
    return genome.split("/")[-1].replace('.fasta','').replace('.faa','').replace('.fna','').replace('.fa','').replace('.FA','').replace('FAA','').replace('.FNA','').replace('.FASTA','').replace('.gff','').replace('.annotated_gff','')