# uniref_annotator.py

This is a script to annotate a fasta file of coding sequences against HUMAnN2-
formatted UniRef90/UniRef50 databases. The script performs a DIAMOND search of 
the input against the databases, and then enforces UniRef's clustering criteria 
on the mapping results to annotate the input sequences. Specifically, to receive 
a UniRef90 annotation, a sequence must have 1) >=90% alignment identify to and
2) >=80% mutual coverage of the centroid.

This approximates UniRef's clustering criteria in which a cluster member must 
have X% identity (e.g. 90 in the case of UniRef90) and 80% coverage of the 
cluster seed. The annotation is approximate in that our comparisons are made 
with cluster representatives (i.e. the sequences in the UniRef files), which are
not necessarily the cluster seed (longest) sequences.

**Author**: Eric Franzosa (eric.franzosa@gmail.com)

## Requirements

* Python 2.7+
* DIAMOND v0.8+
* UniRef90 and UniRef50 databases formatted for DIAMOND
* Mapping from UniRef90 to UniRef50 representatives (recommended)

## Sample command

```
python uniref_annotator.py \
    my_pangenome.ffn \
    --seqtype cds \
    --uniref90db uniref90.dmnd \
    --uniref50db uniref50.dmnd \
    --diamond-options "--threads 8" \
    --transitive-map map_uniref90_uniref50.tsv \
```

This will produce a file `my_pangenome.ffn.annotated` with updated FASTA headers:

```
OLD:
>xx|xxxxx|xxxxxxx
NEW:
>xx|xxxxx|xxxxxxx|UniRef90_ABC|UniRef50_XYZ
NEW, FAILED MAPPING:
>xx|xxxxx|xxxxxxx|UniRef90_unknown|UniRef50_XYZ
```

When running with `--seqtype cds`, the script will translate coding sequences to
protein and then run a `blastp`-like search. I have found this to produce better 
results than simply `blastx`-ing the CDSs vs. UniRef90/50.

In the event that a CDS maps to a UniRef90, the `--transitive-map` option will
look up that UniRef90's corresponding UniRef50 rather than inferring the UniRef50
from mapping results.

## Runtime

Annotating the *S. cerevisiae* pangenome (~6K sequences; ~10 MB) using the above
command required ~20 minutes of real-world time (using 8 cores for search). 
Running the same command on existing mapping results takes ~20 seconds.

## Database paths

As of the time of writing (10/2017), you can find compatible database files at the 
following locations:

```
uniref90.dmnd -> /n/huttenhower_lab/data/humann2_databases/uniref_annotated/uniref90/v1.1_uniref90/uniref90_annotated.1.1.dmnd
uniref50.dmnd -> /n/huttenhower_lab/data/humann2_databases/uniref_annotated/v1.1_uniref50/uniref50_annotated.1.1.dmnd
map_uniref90_uniref50.tsv -> /n/huttenhower_lab/data/idmapping/map_UniRef90_UniRef50.uniq.dat
```

**NOTE**: The above DIAMOND databases are formatted for DIAMOND v0.8, and will not be compatible with DIAMOND v0.9+.