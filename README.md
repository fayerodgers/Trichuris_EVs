# Trichuris_EVs
Scripts for analysis of organoid/EV interactions

## Software and data releases used

kallisto 0.43.1

R version 3.5.0

Sylamer 

Mouse cDNA FASTA file, Ensembl release 97: <ftp://ftp.ensembl.org/pub/release-97/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz>

Mouse gene to transcript mapping, retrieved from Ensembl release 97 BioMart (GRCm38.p6): transcripts2genes.E97.tsv


## Workflow

### Differential expression

experiment_design.tsv contains sample IDs and meta data.

Use the pf scripts to retrieve the FASTQ files from irods, making symlinks in the working directory:

```
tail -n +2 experiment_design.tsv | cut -f 1 | pf data -i - -t file --filetype fastq --file-id-type sample --symlink . --rename
```

Quantify genes counts with Kallisto:

```
tail -n +2 experiment_design.tsv | cut -f 2 | while read id ; do mkdir ${id} ; bsub -o ${id}.o -e ${id}.e -R 'select[mem>=5000] rusage[mem=5000] span[hosts=1]' -M 5000 kallisto quant -i /lustre/scratch118/infgen/team133/fr7/trichuris_projects/Mus_musculus.GRCm38.cdna.E97.all.idx -o ${id} -b 100 ${id}_1.fastq.gz ${id}_2.fastq.gz; done
```

Make a table of reads and Kallisto counts (pseudoaligned reads) - useful for plotting later:

```
tail -n +2 experiment_design.tsv | cut -f 2 | while read id ; do reads=$(zcat ${id}_1.fastq.gz |  awk '{r++}END{print r/4}' )  ; counts=$(awk '{c+=$4}END{print c}' ${id}/abundance.tsv) ;echo ${id}$'\t'${reads}$'\t'${counts} >> counts_table.tsv ; done
```

The R script for differential expression, and generating plots is ```Trichuris_EVs.R```

Genes for plotting in heat maps were selected as follows:

GO term enrichment of regulated genes (padj<0.05) carried out in InnateDB v5.4 (hypergeometric method with BH correction).

Genes associated with the GO terms GO:0051607/GO:0009615 (defense response to virus/response to virus), GO:0070062 (extracellular vesicular exosome) and GO:0045087 (innate immune response) were parsed from the InnateDB output file as follows (eg):

```
cut -f 2,9 InnateDB_goanalysis_0_15_005.txt | egrep "^GO:0051607|^GO:0009615" | cut -f 2 | tr -d '[:space:]' | tr ";" "\n" | sort | uniq > virus.txt 
```

### Motif enrichment in 3' UTRs of differentially expressed genes

Fetch all mouse 3' UTR sequences from Ensembl (using the Perl API) with ```fetch_longest_utrs.pl``` (for each gene, retrieves the 3' UTR sequence of all transcripts and returns the longest).

Filter the file for UTRs of genes that are expressed in our conditions:

```
awk '$2>0{print $1}' deseq2_results.tsv | sed -e 's/"//g' | while read id; do grep -A1 $id longest_utrs.fa >> expressed_utrs.fa ; done
```

Sort the file from shortest UTR to longest:

```
cat expressed_utrs.fa | sed -e '$!N;s/\n/\t/' | while read line ; do echo $(echo $line | cut -f 2 | wc -c)$'\t'$line; done | sort -nk1,1 | cut -f 2,3 | tr '[:space:]' '\n'
```

