# Trichuris_EVs
Scripts for analysis of organoid/EV interactions

## Software and data releases used

kallisto 0.43.1

Mouse cDNA FASTA file, Ensembl release 97: ftp://ftp.ensembl.org/pub/release-97/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz


## Workflow

### Differential expression

experiment_design.tsv contains sample IDs and meta data.

Use the pf scripts to retrieve the FASTQ files from irods, making symlinks in the working directory:

tail -n +2 experiment_design.tsv | cut -f 1 | pf data -i - -t file --filetype fastq --file-id-type sample --symlink . --rename

Quantify genes counts with Kallisto:

tail -n +2 experiment_design.tsv | cut -f 2 | while read id ; do mkdir ${id} ; bsub -o ${id}.o -e ${id}.e -R 'select[mem>=5000] rusage[mem=5000] span[hosts=1]' -M 5000 kallisto quant -i /lustre/scratch118/infgen/team133/fr7/trichuris_projects/GRCm38.cdna.idx -o ${id} -b 100 ${id}_1.fastq.gz ${id}_2.fastq.gz; done


