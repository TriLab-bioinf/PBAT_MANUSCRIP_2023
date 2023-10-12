#!/bin/bash

#############################################
# To run the script
# for SAMPLE in $(cat sample_prefix.txt); do
# ./process_forward_reads.sh ${SAMPLE}
# done
#############################################

SAMPLE=$1

echo
echo "#############################"
echo "#" Trimming ${SAMPLE}
echo "#############################"
echo

# Trim reads
TRIMDIR=cutadapt_trim
if [[ ! -d ./${TRIMDIR} ]]; then
    mkdir ${TRIMDIR}
fi

cutadapt -b "file:primer_A.fasta" \
    --overlap 6 \
    -m 27 \
    --discard-untrimmed ./raw_fastq/${SAMPLE}_L001_R1_001.fastq.gz | \
    cutadapt -b "file:primer_B.fasta" \
        --overlap 6  \
        -m 27 - | \
    cutadapt -b "file:illumina_adapters.fasta" \
        --times 2  \
        --overlap 6  \
        -m 25 \
        -o ./cutadapt_trim/${SAMPLE}_step3.fastq - 1>${SAMPLE}.cutadapt.log 2>&1

# Map reads to genome
BAMDIR=bam
if [[ ! -d ./${BAMDIR} ]]; then
    mkdir ${BAMDIR}
fi

echo
echo "#############################"
echo "#" Mapping ${SAMPLE}
echo "#############################"
echo

bowtie2 --end-to-end \
    --time \
    --threads 8  \
    -x ./GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index \
    -U ./cutadapt_trim/${SAMPLE}_step3.fastq | \
        samtools view -hb - | \
        samtools sort -T ${SAMPLE}.tmp \
            -O BAM \
            -@ 8 \
            --write-index \
            -o ./${BAMDIR}/${SAMPLE}.cutadap_step3.R1.bam - 1>${SAMPLE}.bowtie2.log 2>&1

# Identify peaks with macs2

echo
echo "######################################"
echo "#" Looking for peaks on  sample
echo "#" ${SAMPLE}
echo "######################################"
echo

macs2 callpeak -g 'hs' \
    -t ./${BAMDIR}/${SAMPLE}.cutadap_step3.R1.bam \
    --outdir ./macs2 \
    -n ${SAMPLE} \
    -B \
    --nomodel \
    --extsize 30 \
    --shift -2 1>${SAMPLE}.macs2.log 2>&1


# Convert bam to bed with bedtools
echo
echo "######################################"
echo "#" Converting bam 2 bed
echo "#" ${SAMPLE}
echo "######################################"
echo


BEDDIR=bed

if [[ ! -d ./${BEDDIR} ]]; then
    mkdir ${BEDDIR}
fi

bedtools bamtobed -cigar -i ./${BAMDIR}/${SAMPLE}.cutadap_step3.R1.bam > ./${BEDDIR}/${SAMPLE}.cutadap_step3.R1.bed


# Generate histograms from bed
echo
echo "######################################"
echo "#" Generating histograms
echo "#" ${SAMPLE}
echo "######################################"
echo

HISTDIR=histograms

if [[ ! -d ./${HISTDIR} ]]; then
    mkdir ${HISTDIR}
fi

cat ./${BEDDIR}/${SAMPLE}.cutadap_step3.R1.bed | ./bed2histogram.pl > ./${HISTDIR}/${SAMPLE}.step3.R1.hist



# Extract flanking DNA sequences
echo
echo "######################################"
echo "#" Extracting flanking sequences
echo "#" ${SAMPLE}
echo "######################################"
echo

FLANKDIR=flanking

if [[ ! -d ./${FLANKDIR} ]]; then
    mkdir ${FLANKDIR}
fi

seqtk subseq human_GRCh38.fasta ./histograms/${SAMPLE}.step3.R1.hist > ./${FLANKDIR}/${SAMPLE}.IS_flank.fasta


echo
echo "######################################"
echo "#" ${SAMPLE} Done !
echo "######################################"
echo ------------------------------------------------------------------------------------------------
echo