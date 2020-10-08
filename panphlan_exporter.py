
__author__ = ('Leonard Dubois (leonard.dubois@unitn.it), '
              'Aitor Blanco (aitor.blancomiguez@unitn.it)')
__version__ = '0.01'
__date__ = '6 Oct 2020'


import os, bz2, subprocess, sys, argparse, time, shutil, gffutils
from glob import iglob
from BCBio import GFF
from Bio import SeqIO
from utils import info, error, get_genome_name
from parallelisation import execute_pool
from external_exec import execute_bowtie2_build, execute_bowtie2_inspect, execute_prokka, execute_uniref_annotator, execute_mmseq2
from external_exec import UNIREF_VERSION

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

    p.add_argument('-d', '--db_path', type=str, default='.',
                   help="Path to UniRef DIAMOND formatted database")

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



# ------------------------------------------------------------------------------
#            Checking softwares
# ------------------------------------------------------------------------------

def check_software(software_name):
    platform = sys.platform.lower()[0:3]

    if platform == 'win':
        p = subprocess.Popen(['where', software_name], stdout=subprocess.PIPE)
    else: # Linux, Mac, ...
        p = subprocess.Popen(['which', software_name], stdout=subprocess.PIPE)
    p_out = p.communicate()[0].decode('utf-8')
    if not p_out == '':
        print('[I] {} is installed'.format(software_name))
        print('   path: {}'.format(p_out.strip()))
        return(p_out.strip())
    else:
        print('\n[E] {} software not found.\n'.format(software_name))
        print('\n[E] Please, install {} or update your PATH variable.\n'.format(software_name))
        sys.exit()



# ------------------------------------------------------------------------------
#            Fuctions retrieved from generate pangenome script
# ------------------------------------------------------------------------------

def parallel_prokka(genome, outputdir):
    genome_name = get_genome_name(genome)
    os.mkdir(os.path.join(outputdir, 'prokka', genome_name))
    execute_prokka(genome, os.path.join(outputdir, 'prokka', genome_name), genome_name)

def parallel_uniref_annotator(genome, outputdir, db_path, diamond_path):
    genome_name = get_genome_name(genome)
    os.mkdir(os.path.join(outputdir, 'uniref', genome_name))
    os.mkdir(os.path.join(outputdir, 'uniref', genome_name, 'tmp'))
    input_file = os.path.join(outputdir, 'prokka', genome_name, genome_name + '.faa')
    execute_uniref_annotator(input_file, os.path.join(outputdir, 'uniref', genome_name, genome_name + '.faa'), os.path.join(outputdir, 'uniref', genome_name, 'tmp'), db_path, diamond_path)
    shutil.rmtree(os.path.join(outputdir, 'uniref', genome_name, 'tmp'))


# CLUSTER UNIREF -----------------------------------------------------------

def get_clusters_90(outputdir):
    clusters = list()
    with open(os.path.join(outputdir, 'unannotated_90.clustered.tsv'), 'r') as tsv_reader:
        cluster_ref = ''
        for row in tsv_reader:
            row = row.strip().split('\t')
            if cluster_ref == row[0] and row[0] != row[1]:
                cluster += [row[1]]
            else:
                if cluster_ref != '':
                    clusters.append({'cluster': 'UniRef90_UNK-' + cluster_ref, 'proteins': cluster})
                cluster_ref = row[0]
                cluster = [row[0]]
                if row[0] != row[1]:
                    cluster += [row[1]]
    clusters.append({'cluster': 'UniRef90_UNK-' + cluster_ref, 'proteins':cluster})
    return clusters


def get_fasta_uniref_level(uniref_level):
    if uniref_level == '50':
        return -1
    elif uniref_level == '90':
        return -2
    else:
        raise ValueError('UniRef level not supported.')


def get_unannotated_proteins(inputdir, outputdir, uniref_level):
    if uniref_level == '50':
        uniref90_clusters = list()
        for cluster in get_clusters_90(outputdir):
            uniref90_clusters.append(cluster['cluster'][13:])

    with open(os.path.join(outputdir, 'unannotated', 'unannotated_' + uniref_level +'.faa'),'w') as unnanotated_proteins:
        for genome in iglob(inputdir + '/*'):
            genome_name = get_genome_name(genome)
            write_fasta = False
            with open(os.path.join(outputdir, 'uniref', genome_name, genome_name + '.faa'), 'r') as read_fasta:
                for line in read_fasta:
                    if line.startswith('>'):
                        if line.strip().split('|')[get_fasta_uniref_level(uniref_level)].split('_')[1] == 'unknown':
                            if uniref_level == '90' or (uniref_level == '50' and line.strip().split(' ')[0][1:] in uniref90_clusters):
                                write_fasta = True
                            else:
                                write_fasta = False
                            if write_fasta:
                                unnanotated_proteins.write(line)
                        else:
                            write_fasta = False
                    elif write_fasta:
                        unnanotated_proteins.write(line)


def create_mmseqs_db(input_fasta, output_dir, mmseq2_path):
    file_name = os.path.splitext(os.path.basename(input_fasta))[0]
    execute_mmseq2('createdb '+input_fasta+' '+os.path.join(output_dir, 'db', file_name))


def get_uniclust_params(uniref_level):
    if uniref_level == '50':
        return ('0.8', '0.5')
    elif uniref_level == '90':
        return ('0.8', '0.9')
    else:
        raise ValueError('UniRef level not supported.')

def cluster_db(input_fasta, output_dir, uniref_level, nprocs, mmseq2_path):
    file_name = os.path.splitext(os.path.basename(input_fasta))[0]
    coverage, identity = get_uniclust_params(uniref_level)
    execute_mmseq2('cluster ' + os.path.join(output_dir, 'db', file_name) + ' ' +
        os.path.join(output_dir, 'db_clustered', file_name) + ' ' +
        os.path.join(output_dir,'tmp') + ' -c ' + coverage + ' --min-seq-id ' + identity +
        ' --threads ' + nprocs)

def create_mmseqs_folders(outputdir):
    os.mkdir(os.path.join(outputdir, 'mmseq'))
    os.mkdir(os.path.join(outputdir, 'mmseq', 'tmp'))
    os.mkdir(os.path.join(outputdir, 'mmseq', 'db'))
    os.mkdir(os.path.join(outputdir, 'mmseq', 'db_clustered'))


def get_clustering_results(input_fasta, output_dir, uniref_level, nprocs):
    file_name = os.path.splitext(os.path.basename(input_fasta))[0]
    execute_mmseq2('createtsv ' + os.path.join(output_dir, 'mmseq', 'db', file_name) + ' ' +
        os.path.join(output_dir, 'mmseq', 'db', file_name) + ' ' +
        os.path.join(output_dir, 'mmseq', 'db_clustered', file_name) + ' ' +
        os.path.join(output_dir, file_name + '.clustered.tsv') + ' --threads ' + nprocs)
    os.remove(input_fasta)
    return os.path.join(output_dir, file_name + '_clustered_' + uniref_level + '.tsv')


def cluster_unnanotated_proteins(input_fasta, output_dir, uniref_level, nprocs, mmseq2_path):
    create_mmseqs_folders(output_dir)
    create_mmseqs_db(input_fasta, os.path.join(output_dir, 'mmseq'), mmseq2_path)
    cluster_db(input_fasta, os.path.join(output_dir, 'mmseq'), uniref_level, str(nprocs), mmseq2_path)
    clustering_file = get_clustering_results(input_fasta,  output_dir, uniref_level, str(nprocs))
    shutil.rmtree(os.path.join(output_dir, 'mmseq'))
    return clustering_file


# RE ANNOTATE GENOME --------------------------------

def iterator_gff(db, nr_faa_fp):
    with open(nr_faa_fp) as nr_aa:
        for entry in iter(entry.strip() for entry in nr_aa if entry.startswith('>')):
            cds_id = entry.split(' ')[0][1:]
            nr90, nr50 = entry.split('|')[1:]
            gene = db[cds_id]
            attrs = dict(gene.attributes)
            attrs.update({'UniRef90' : [nr90], 'UniRef50' : [nr50]})
            attrs = gffutils.attributes.Attributes(attrs)
            gene.attributes = attrs

            yield gene

def update_gff_nr_annot(genome, outputdir):
    genome_name = get_genome_name(genome)
    prokka_gff_fp = os.path.join(outputdir, 'prokka', genome_name, genome_name + '.gff')
    nr_faa_fp = os.path.join(outputdir, 'uniref', genome_name, genome_name + '.faa')

    prokka_gff = gffutils.create_db(prokka_gff_fp, ':memory:', id_spec='ID', merge_strategy='merge')
    results_update = iterator_gff(prokka_gff, nr_faa_fp)
    prokka_gff.update(results_update, merge_strategy='replace')

    with open(prokka_gff_fp.replace('.gff', '.annotated_gff'), 'w') as fout, \
            open(prokka_gff_fp) as fna_in:
        fout.write('\n'.join('##{}'.format(d) for d in prokka_gff.directives))
        fout.write('\n')
        fout.write('\n'.join(str(feature) for feature in prokka_gff.all_features()))
        fout.write('\n##FASTA\n')
        fout.write(''.join(contig.format('fasta') for contig in SeqIO.parse(fna_in, 'fasta')))

    [ prokka_gff.delete(x) for x in prokka_gff.all_features() ]


def get_uniref90_to_50(inputdir, outputdir, DATABASES_PATH):
    uniref90_to_50 = dict()
    k90_u50 = list()
    for genome in iglob(inputdir + '/*'):
        genome_name = get_genome_name(genome)
        with open(os.path.join(outputdir, 'prokka', genome_name, genome_name + '.annotated_gff'), 'r') as gff_file:
            for rec in gff_file:
                if 'UniRef50_unknown' in rec and not 'UniRef90_unknown' in rec:
                    k90_u50.append('UniRef90_' + rec.split('UniRef90_')[1].split(';')[0])
    if len(k90_u50) > 0:
        with open('{}/UniRef90to50_{}.tsv'.format(DATABASES_PATH, UNIREF_VERSION),'r') as read_file:
            for line in read_file:
                line = line.strip().split('\t')
                if line[0] in k90_u50:
                    uniref90_to_50[line[0]] = line[1]
                    k90_u50.remove(line[0])
                    if len(k90_u50) == 0:
                        break
    return uniref90_to_50


def get_clusters_50(inputdir, outputdir, clusters_90):
    uniref90_clusters = dict()
    for cluster in clusters_90:
        uniref90_clusters[cluster['cluster'][13:]] = cluster['proteins']

    clusters = list()
    with open(os.path.join(outputdir, 'unannotated_50.clustered.tsv'), 'r') as tsv_reader:
        cluster_ref = ''
        for row in tsv_reader:
            row = row.strip().split('\t')
            if cluster_ref == row[0] and row[0] != row[1]:
                cluster += uniref90_clusters[row[1]]
            else:
                if cluster_ref != '':
                    clusters.append({'cluster': 'UniRef50_UNK-' + cluster_ref, 'proteins': cluster})
                cluster_ref = row[0]
                cluster = uniref90_clusters[row[0]]
                if row[0] != row[1]:
                    cluster += uniref90_clusters[row[1]]
    clusters.append({'cluster': 'UniRef50_UNK-' + cluster_ref, 'proteins':cluster})
    return clusters


def get_u90_to_k50(inputdir, outputdir):
    u90_to_k50 = dict()
    for genome in iglob(inputdir + '/*'):
        genome_name = get_genome_name(genome)
        with open(os.path.join(outputdir, 'uniref', genome_name, genome_name + '.faa'), 'r') as read_fasta:
            for line in read_fasta:
                if line.startswith('>'):
                    if line.strip().split('|')[-2].split('_')[1] == 'unknown' and line.strip().split('|')[-1].split('_')[1] != 'unknown':
                        u90_to_k50[line.split(' ')[0][1:]] = line.strip().split('|')[-1]
    return u90_to_k50



def get_clusters_by_file(clusters, prokka_id):
    clusters_by_file = list()
    for c in clusters:
        proteins_by_cluster = []
        for p in c['proteins']:
            if p.split('|')[0].startswith(prokka_id):
                proteins_by_cluster.append(p.split('|')[-1])
        if len(proteins_by_cluster) > 0:
            clusters_by_file.append({'cluster': c['cluster'],
                'proteins': proteins_by_cluster})
    return clusters_by_file



def get_prokka_id(genome_name, outputdir):
    with open(os.path.join(outputdir, 'prokka', genome_name, genome_name + '.faa'), 'r') as read_fasta:
        return read_fasta.readline().split(' ')[0].split('_')[0][1:]


def cluster_by_file(clusters_50, clusters_90, outputdir):
    clustering_by_file = []
    for f in iglob(os.path.join(outputdir, 'prokka/*/*.annotated_gff')):
        genome_name = get_genome_name(f)
        prokka_id = get_prokka_id(genome_name, outputdir)
        clusters_50_by_file = get_clusters_by_file(clusters_50, prokka_id)
        clusters_90_by_file = get_clusters_by_file(clusters_90, prokka_id)
        clustering_by_file.append({'gff': f, 'clusters_50': clusters_50_by_file, 'clusters_90': clusters_90_by_file})
    return clustering_by_file


def fix_u90_to_k50_inconsistencies(gff_list, u90_to_k50):
    for n, rec in enumerate(gff_list):
        if 'UniRef50_unknown' in rec and 'UniRef90_UNK-' in rec:
            u90_id = rec.split('UniRef90_UNK-')[1].split(';')[0]
            if u90_id in u90_to_k50:
                gff_list[n] = rec.replace('UniRef50_unknown', u90_to_k50[u90_id])
    return gff_list


def get_reannotated_gff(gff_list, clusters, uniref90_to_50, uniref_level):
    if len(clusters) > 0:
        for n, rec in enumerate(gff_list):
            if uniref_level == '50' and 'UniRef50_unknown' in rec and not 'UniRef90_unknown' in rec:
                gff_list[n] = rec.replace('UniRef50_unknown', uniref90_to_50['UniRef90_' + rec.split('UniRef90_')[1].split(';')[0]])
            else:
                for c in clusters:
                    for p in c['proteins']:
                        if 'UniRef' + uniref_level + '_unknown' in rec and p == rec.split('\t')[8].split(';')[0].replace('ID=',''):
                            gff_list[n] = rec.replace('UniRef' + uniref_level + '_unknown', c['cluster'])
                            c['proteins'].remove(p)
                            break
    return gff_list


def reannotate_gff(file_clustered, output_dir, uniref90_to_50, u90_to_k50):
    initial_gff = file_clustered['gff']
    clusters_50 = file_clustered['clusters_50']
    clusters_90 = file_clustered['clusters_90']
    gff_list = list()
    with open(initial_gff, 'r') as gff_file:
        gff_list = [line for line in gff_file]
        with bz2.open(os.path.join(output_dir, 'annotations', get_genome_name(initial_gff) + '.gff.bz2'), 'wt') as output_gff:
            gff_list = get_reannotated_gff(gff_list, clusters_50, uniref90_to_50, '50')
            gff_list = get_reannotated_gff(gff_list, clusters_90, uniref90_to_50, '90')
            gff_list = fix_u90_to_k50_inconsistencies(gff_list, u90_to_k50)
            for rec in gff_list:
                output_gff.write(rec)

def reannotate_genomes(inputdir, outputdir, db_path):
    for genome in iglob(inputdir + '/*'):
        update_gff_nr_annot(genome, outputdir)
    uniref90_to_50 = get_uniref90_to_50(inputdir, outputdir, db_path)
    clusters_90 = get_clusters_90(outputdir)
    clusters_50 = get_clusters_50(inputdir, outputdir, clusters_90)
    u90_to_k50 = get_u90_to_k50(inputdir, outputdir)
    clustering_by_file = cluster_by_file(clusters_50, clusters_90, outputdir)
    for file_clustered in clustering_by_file:
        reannotate_gff(file_clustered, outputdir, uniref90_to_50, u90_to_k50)



# ------------------------------------------------------------------------------
#           PanPhlAn exporting functions
# ------------------------------------------------------------------------------

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


def create_bt2_indexes(ppa_outdir, clade_name, nprocs, bowtie2_path):
    fna_file = os.path.join(ppa_outdir, clade_name + '_contigs.fna')
    indexes_path = os.path.join(ppa_outdir, clade_name)
    execute_bowtie2_build(fna_file, indexes_path, nprocs, bowtie2_path)
    execute_bowtie2_inspect(indexes_path, bowtie2_path)


def create_temporal_folders(tmp_dir):
    os.mkdir(os.path.join(tmp_dir, 'tmp'))
    os.mkdir(os.path.join(tmp_dir, 'tmp', 'prokka'))
    os.mkdir(os.path.join(tmp_dir, 'tmp', 'uniref'))
    os.mkdir(os.path.join(tmp_dir, 'tmp', 'annotations'))
    os.mkdir(os.path.join(tmp_dir, 'tmp', 'unannotated'))
    return os.path.join(tmp_dir, 'tmp')


# ------------------------------------------------------------------------------
#                    MAIN FUNCTION
# ------------------------------------------------------------------------------

def panphlan_exporter(inputdir, tmp_dir, ppa_outdir, clade_name, nprocs, db_path):
    mmseq2_path = check_software("mmseqs")
    diamond_path = check_software("diamond")
    bowtie2_path = check_software("bowtie2")
    extend_pangenome = extend_panphlan_output(ppa_outdir)
    if extend_pangenome:
        info("PanPhlAn pangenome already exists. It will be extended with genomes provided")
    info("Creating temporal folders...")
    tmp_dir = create_temporal_folders(tmp_dir)
    info("Done.")
    info("Checking PanPhlAn output directory...")
    contigs_names_dict = get_contig_names_dict(inputdir)
    info("Done.")
    info("Executing prokka...")
    execute_pool(((parallel_prokka, genome, tmp_dir) for genome in iglob(inputdir + "/*")), max(1, int(args.nprocs/4)))
    info("Done.")
    info("Executing uniref_annotator...")
    execute_pool(((parallel_uniref_annotator, genome, tmp_dir, db_path, diamond_path) for genome in iglob(inputdir + "/*")), max(1, int(args.nprocs/4)))
    info("Done.")
    info("Clustering unnanotated proteins at UniRef90 level...")
    get_unannotated_proteins(inputdir, tmp_dir, '90')
    cluster_unnanotated_proteins(os.path.join(tmp_dir, 'unannotated', 'unannotated_90.faa'), tmp_dir, '90', nprocs, mmseq2_path)
    info('Done.')
    info('Clustering unnanotated proteins at UniRef50 level...')
    get_unannotated_proteins(inputdir, tmp_dir, '50')
    cluster_unnanotated_proteins(os.path.join(tmp_dir, 'unannotated', 'unannotated_50.faa'), tmp_dir, '50', nprocs, mmseq2_path)
    info('Done.')
    info('Reannotating genomes...')
    reannotate_genomes(inputdir, tmp_dir, db_path)
    contigs_names_dict_prokka = get_contig_names_dict_prokka(inputdir, tmp_dir)
    info('Done.')
    info("Writing PanPhlAn tsv...")
    write_panphlan_tsv(inputdir, tmp_dir, ppa_outdir, clade_name, contigs_names_dict, contigs_names_dict_prokka, extend_pangenome)
    info("Done")
    info("Writing PanPhlAn fna...")
    write_panphlan_contigs(inputdir, ppa_outdir, clade_name )
    info("Done")
    info("Generating indexes for bowtie2...")
    create_bt2_indexes(ppa_outdir, clade_name, nprocs, bowtie2_path)
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
    panphlan_exporter(args.input, args.tmp, args.output, args.clade_name, args.nprocs, args.db_path)

    exec_time = time.time() - t0
    info('Finish execution ({} seconds)\n'.format(round(exec_time, 2)))
