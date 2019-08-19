# Trichuris_EVs
Scripts for analysis of organoid/EV interactions

## Software and data releases used

kallisto 0.43.1

R 3.5.0

dustmasker 1.0.0

mkvtree 2.3.0, vmatch 2.3.0 (run from purge-sequence in RSAT)

Sylamer 18-131

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
cat expressed_utrs.fa | sed -e '$!N;s/\n/\t/' | while read line ; do echo $(echo $line | cut -f 2 | wc -c)$'\t'$line; done | sort -nk1,1 | cut -f 2,3 | tr '[:space:]' '\n' > expressed_utrs_sorted.fa
```

Dust and purge the UTRs (then hard mask):

```
dustmasker -in expressed_utrs_sorted.fa -out expressed_utrs_dusted.fa -outfmt fasta
purge-sequence -i expressed_utrs_dusted.fa  -1str -ml 40 -mis 2 -skip_short 20 -o expressed_utrs_dusted_purged.fa 
sed -e 's/[a-z]/N/g' expressed_utrs_dusted_purged.fa > expressed_utrs_hard_masked.fa 
```

Generate a rankfile of Trichuris genes, from most significantly down-regulated to most significantly up-regualated:

```
tail -n +2 deseq2_results.tsv | sed -e 's/NA/0.99/g' | sed -e 's/"//g'| awk '$2>0 {print $1,$3,-log($6)/log(10)}'| while read -r id change logpval; do if [[ "$change" =~ '-' ]]; then logpval='-'$logpval ; fi ; echo $id$'\t'$logpval ; done | sort -g -k2,2 > rankfile.tsv
```

Also generate a rankfile of only the significantly (padj<0.05) downegulated genes:

```
tail -n +2 deseq2_results.tsv | sed -e 's/"//g'| awk '$3<0 && $7<0.05 {print $1,log($6)/log(10)}'| sort -g -k2,2 > subset.tsv
```

Run Sylamer:

```
sylamer -fasta expressed_utrs_hard_masked.fa -universe rankfile.tsv -k 6 -m 4 -grow 50 -o out.universe.6mer
sylamer -fasta expressed_utrs_hard_masked.fa -universe rankfile.tsv -k 7 -m 4 -grow 50 -o out.universe.7mer
sylamer -fasta expressed_utrs_hard_masked.fa -universe rankfile.tsv -k 8 -m 4 -grow 50 -o out.universe.8mer

sylamer -fasta expressed_utrs_hard_masked.fa -subset subset.tsv -k 6 -m 4 -grow 10 -o out.subset.6mer
sylamer -fasta expressed_utrs_hard_masked.fa -subset subset.tsv -k 7 -m 4 -grow 10 -o out.subset.7mer
sylamer -fasta expressed_utrs_hard_masked.fa -subset subset.tsv -k 8 -m 4 -grow 10 -o out.subset.8mer
```

Plot most significant words (eg for k=6):

```
 R --vanilla --slave --quiet --args --data=out.universe.6mer --title='Universe.k6' --pdf=6mer.universe.pdf < ~/bio_software/sylamer/wrapperscript/syl.plot.R
```

Plot selected motifs (eg for k=6):
```
R --vanilla --slave --quiet --args --data=out.universe.6mer --add="ACTGGA,ACTCGT,GACCTG,TTCGAG,CCGTTC,GCTTAG,CTGATC,CTCCAT" --top=0 --bottom=0 --title='Universe.k6' --pdf=6mer.universe.selected.pdf < ~/bio_software/sylamer/wrapperscript/syl.plot.R
[1] "ACTGGA" "ACTCGT" "GACCTG" "TTCGAG" "CCGTTC" "GCTTAG" "CTGATC" "CTCCAT"
```







