
__author__ = ('Leonard Dubois (leonard.dubois@unitn.it), '
              'Aitor Blanco (aitor.blancomiguez@unitn.it)')
__version__ = '0.01'
__date__ = '20 Aug 2020'


import os, bz2, subprocess, sys, argparse, time, shutil
from glob import iglob
from BCBio import GFF
from Bio import SeqIO
from utils import info, error, get_genome_name
from parallelisation import execute_pool
from external_exec import execute_bowtie2_build, execute_bowtie2_inspect
from generate_pangenome import parallel_prokka, parallel_uniref_annotator, get_unannotated_proteins, cluster_unnanotated_proteins
from generate_pangenome import reannotate_genomes

def read_params():
    p = argparse.ArgumentParser(description='')
    p.add_argument('-i', '--input', type=str, default=None,
                   help='The input directory containing the FASTA files')

    p.add_argument('-t', '--tmp', type=str, default=None,
                   help='The temporal output directory. Default: output directory')
                   
    p.add_argument('-p', '--output', type=str, default=None,
                   help="The output directory for PanPhlAn pangenome")
                   
    p.add_argument('-c', '--clade_name', type=str, default='panplhan_clade',
                   help="Prefix of the PanPhlAn output files")
                   
    p.add_argument('-n', '--nprocs', type=int, default=1,
                   help='The number of threads to use')
    
    return p.parse_args()


def check_params(args):
    if not args.input:
        error('-i (or --input) must be specified', exit=True, 
            init_new_line=True)
    elif not args.output:
        error('-o (or --output) must be specified', exit=True, 
            init_new_line=True)
    elif args.tmp is None:
        args.tmp = args.output
    if not os.path.exists(args.output):
        os.mkdir(args.output)    

# PanPhlAn exporting functions -------------------------------------------------


def extend_panphlan_output(ppa_outdir):
    if not os.path.exists(ppa_outdir):
        os.mkdir(ppa_outdir)
        return False
    elif len(os.listdir(ppa_outdir)) >= 8:
        return True
    else:
        return False
        
        
def get_contig_names_dict(inputdir):
    contigs_names_dict = {}
    for f in os.listdir(inputdir):
        genome = get_genome_name(f)
        with open(os.path.join(inputdir, f), 'r') as IN:
            records = list(SeqIO.parse(IN, 'fasta'))
        contigs_names_dict[genome] = []
        for r in records:
            contigs_names_dict[genome].append(r.id)
    return contigs_names_dict


def get_contig_names_dict_prokka(inputdir, tmp_dir):
    contigs_names_dict_prokka = {}
    for genome in iglob(inputdir + '/*'):
        genome_name = get_genome_name(genome)
        contigs_names_dict_prokka[genome_name] = []
        with bz2.open(os.path.join(tmp_dir, 'annotations', genome_name + '.gff.bz2'), 'rt') as IN:
            IN.readline()
            for l in IN:
                if l.startswith('##'):
                    contig_name = l.strip().split(' ')[1]
                    contigs_names_dict_prokka[genome_name].append(contig_name)
                else:
                    break
    return contigs_names_dict_prokka


def convert_contig_name(contig_name_prokka, contigs_names_list, contigs_names_list_prokka):
    idx = contigs_names_list_prokka.index(contig_name_prokka)
    return(contigs_names_list[idx])
    

def write_panphlan_tsv(inputdir, tmp_dir, ppa_outdir, clade_name, contigs_names_dict, contigs_names_dict_prokka, extend_pangenome):
    mode = 'a' if extend_pangenome else 'w'
    OUT_TSV = open(os.path.join(ppa_outdir, clade_name + '_pangenome.tsv'), mode)
    for genome in iglob(inputdir + '/*'):
        genome_name = get_genome_name(genome)
        genome_id = genome_name.split('.')[0]
        gff_file = bz2.open(os.path.join(tmp_dir, 'annotations', genome_name + '.gff.bz2'), 'rt')
        for rec in GFF.parse(gff_file, limit_info=dict(gff_type = ['CDS'])):
            contig_name_prokka = rec.id
            contig_id = convert_contig_name(contig_name_prokka, contigs_names_dict[genome_name], contigs_names_dict_prokka[genome_name])
            for f in rec.features:
                if 'gene' in f.qualifiers:
                    name = f.qualifiers['gene'][0]
                else:
                    name = f.id
                feature = str(f.location)
                init = int(feature.split(':')[0].split('[')[1])
                end = int(feature.split(':')[1].split(']')[0])
                OUT_TSV.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(f.qualifiers['UniRef90'][0], name, genome_id, contig_id, init, end))
    OUT_TSV.close()

    
def write_panphlan_contigs(inputdir, ppa_outdir, clade_name ):
    OUT = open(os.path.join(ppa_outdir, clade_name + '_contigs.fna'), 'a')
    for genome in iglob(inputdir + '/*'):
        with open(genome, 'r') as IN:
            for l in IN:
                OUT.write(l)
    OUT.close()
    

def create_bt2_indexes(ppa_outdir, clade_name, nprocs):    
    fna_file = os.path.join(ppa_outdir, clade_name + '_contigs.fna')
    indexes_path = os.path.join(ppa_outdir, clade_name)
    execute_bowtie2_build(fna_file, indexes_path, nprocs)
    execute_bowtie2_inspect(indexes_path)
    

def create_temporal_folders(tmp_dir):
    os.mkdir(os.path.join(tmp_dir, 'tmp'))
    os.mkdir(os.path.join(tmp_dir, 'tmp', 'prokka'))
    os.mkdir(os.path.join(tmp_dir, 'tmp', 'uniref'))
    os.mkdir(os.path.join(tmp_dir, 'tmp', 'annotations'))
    os.mkdir(os.path.join(tmp_dir, 'tmp', 'unannotated'))
    return os.path.join(tmp_dir, 'tmp')

    
def panphlan_exporter(inputdir, tmp_dir, ppa_outdir, clade_name, nprocs):
    info("Creating temporal folders...")
    tmp_dir = create_temporal_folders(tmp_dir)
    info("Done.")
    info("Checking PanPhlAn output directory...")
    extend_pangenome = extend_panphlan_output(ppa_outdir)
    if extend_pangenome:
        info("PanPhlAn pangenome already exists. It will be extended with genomes provided")
    contigs_names_dict = get_contig_names_dict(inputdir)
    info("Done.")
    info("Executing prokka...")
    execute_pool(((parallel_prokka, genome, tmp_dir) for genome in iglob(inputdir + "/*")), max(1, int(args.nprocs/4)))
    info("Done.")
    info("Executing uniref_annotator...")
    execute_pool(((parallel_uniref_annotator, genome, tmp_dir) for genome in iglob(inputdir + "/*")), max(1, int(args.nprocs/4)))
    info("Done.")
    info("Clustering unnanotated proteins at UniRef90 level...")
    get_unannotated_proteins(inputdir, tmp_dir, '90')
    cluster_unnanotated_proteins(os.path.join(tmp_dir, 'unannotated', 'unannotated_90.faa'), tmp_dir, '90', nprocs)
    info('Done.')
    info('Clustering unnanotated proteins at UniRef50 level...')
    get_unannotated_proteins(inputdir, tmp_dir, '50')
    cluster_unnanotated_proteins(os.path.join(tmp_dir, 'unannotated', 'unannotated_50.faa'), tmp_dir, '50', nprocs)
    info('Done.')
    info('Reannotating genomes...')
    reannotate_genomes(inputdir, tmp_dir)
    contigs_names_dict_prokka = get_contig_names_dict_prokka(inputdir, tmp_dir)
    info('Done.')
    info("Writing PanPhlAn tsv...")
    write_panphlan_tsv(inputdir, tmp_dir, ppa_outdir, clade_name, contigs_names_dict, contigs_names_dict_prokka, extend_pangenome)
    info("Done")
    info("Writing PanPhlAn fna...")
    write_panphlan_contigs(inputdir, ppa_outdir, clade_name )
    info("Done")
    info("Generating indexes for bowtie2...")
    create_bt2_indexes(ppa_outdir, clade_name, nprocs)
    info("Done")
    info('Removing temporal files...')
    shutil.rmtree(tmp_dir)
    info('Done.')


if __name__ == '__main__':    
    info('Start execution', init_new_line=False)
    t0 = time.time()
    args = read_params()
    check_params(args)
    print(args)
    panphlan_exporter(args.input, args.tmp, args.output, args.clade_name, args.nprocs)

    exec_time = time.time() - t0
    info('Finish execution ({} seconds)\n'.format(round(exec_time, 2)))

# ------------------------------------------------------------------------------