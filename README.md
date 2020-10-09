## UniRef genes families-level pangenome building and annotation

This tools provides a pipeline for annotating and clustering input genomes sequences into UniRef90/UniRef50 genes families and clustering unknown coding sequences.
The output provided is a ready-to-use PanPhlAn pangenome. Thus, it will countain all genomes contigs in a multi-FASTA file, precomputed bowtie2 indexes, and a pangenome tsv file mapping gene location on contigs.

### Pipeline



### Dependencies :

The following Python packages are needed .
* BioPython
* bcbio-gff
* gffutils

The following external tools should be installed (and the PATH variable properly configured) :
* Prokka (https://github.com/tseemann/prokka)
* MMSEQ2 (https://github.com/soedinglab/MMseqs2)
* DIAMOND (https://github.com/bbuchfink/diamond)
* BowTie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

On top on that, UniRef DIAMOND databases should be downloaded via the `download_databases.py` script.


### Usage

```
python panphlan_exporter.py --input [input_genomes_folder]          \
                            --output [output_pangenome_folder]      \
                            --db_path [path_to_UniRef_DIAMOND_databases]
```
* The `--input  [input_genomes_folder]` should contain one fasta file per genome. The script assumes that the file name is the genome name
* The `--output [output_pangenome_folder]` will be created if not existing

Additionnal parameters could be provided :
* `-t` or `--tmp` specifies another directory for temporary files. Default is the output folder
* `-c` or `--clade_name` specifies a prefix for PanPhlAn output files. The best would be the full species name (e.g. `Escherichia_coli`). Default is `panplhan_clade`
* `-n` or `--nprocs` the number of threads to use.





