 
 
 #!/bin/zsh
 
 GENOME_FASTA=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
 GENOME_ANNOTATION=GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf
 
# Export gene coords from genome annotation gtf file as a bed file
echo "#################################"
echo Exporting gene annotations
echo "#################################"

egrep  --color "\tgene\t" ${GENOME_ANNOTATION} | \
    cut -f 1,4,5 | \
    perl -lane '$F[1]--;print "$F[0]\t$F[1]\t$F[2]"' | \
    egrep "^chr\d+\t" > gene_unmerged.bed

# Merge overlapping exons
bedtools merge -i gene_unmerged.bed > gene_merged.bed

# Generate bed file with intergenic regions
echo "#################################"
echo Exporting intergenic annotations
echo "#################################"

seqtk comp ${GENOME_FASTA} | \
    egrep "^chr\d+\t" | \
    perl -lane 'print "$F[0]\t0\t$F[1]"' > chromosomes.bed

cut -f 1,3 chromosomes.bed > my.genome

bedtools complement -i  gene_merged.bed -g my.genome > intergenic.bed

# Centromeric coords were extracted from UCSC Table browser
# https://genome.ucsc.edu/cgi-bin/hgTables
# on human assembly Dec 2012 (GRCh38/hg38)
# Group = "Mapping and Sequencing"
# Track = Centromeres 
# Output format = BED
# Output file centromeres.bed

echo "#################################"
echo Sorting centromeres.bed
echo "#################################"

bedtools sort -i centromeres.bed > centromeres.bed

# Remove centromeric regions from intergenic regions
bedtools subtract -a intergenic.bed -b centromeres.bed > intergenic_no_centromere.bed

# Remove centromeric regions from genes
bedtools subtract -a gene_merged.bed -b centromeres.bed > gene_merged_no_centromere.bed

# Export exon coords from genome annotation gtf file as a bed file
echo "#################################"
echo Exporting exonic annotations
echo "#################################"

egrep  --color "\texon\t" ${GENOME_ANNOTATION} | \
    cut -f 1,4,5 | \
    perl -lane '$F[1]--;print "$F[0]\t$F[1]\t$F[2]"' | \
    egrep "^chr\d+\t" | \
    sort -k1,1 -k2,2n > exon_unmerged.bed

# Merge overlaping exons
bedtools merge -i exon_unmerged.bed > exon_merged.bed

# Export CDS coords from genome annotation gtf file as a bed file
echo "#################################"
echo Exporting CDS annotations
echo "#################################"

egrep  --color "\tCDS\t" ${GENOME_ANNOTATION} | \
    cut -f 1,4,5 | \
    perl -lane '$F[1]--;print "$F[0]\t$F[1]\t$F[2]"' | \
    egrep "^chr\d+\t" | \
    sort -k1,1 -k2,2n > cds_unmerged.bed

# Merge overlaping CDSs
bedtools merge -i cds_unmerged.bed > cds_merged.bed

# Generate bed file with intron coords
bedtools subtract -a gene_merged.bed -b exon_merged.bed > introns_unmerged.bed

# Merge overlaping introns
bedtools merge -i introns_unmerged.bed > introns_merged.bed

# Generate bed file with UTR coords
bedtools subtract -a exon_merged.bed -b cds_merged.bed > utr_unmerged.bed

# Merge overlaping UTRs
bedtools merge -i utr_unmerged.bed > utr_merged.bed

# Remove centromeric regions from CDS regions
bedtools subtract -a cds_merged.bed -b centromeres.bed > cds_merged_no_centromere.bed

# Remove centromeric regions from intronic regions
bedtools subtract -a introns_merged.bed -b centromeres.bed > introns_merged_no_centromere.bed

# Remove centromeric regions from UTR regions
bedtools subtract -a utr_merged.bed -b centromeres.bed > utr_merged_no_centromere.bed


# Count total length for each annotation type with and without regions overlaping centromeres
intron_counts=$(./count_total_length_from_bed.pl < introns_merged.bed)
intron_counts_no_centromere=$(./count_total_length_from_bed.pl < introns_merged_no_centromere.bed)

intergenic_counts=$(./count_total_length_from_bed.pl < intergenic.bed)
intergenic_counts_no_centromere=$(./count_total_length_from_bed.pl < intergenic_no_centromere.bed)

cds_counts=$(./count_total_length_from_bed.pl < cds_merged.bed)
cds_counts_no_centromere=$(./count_total_length_from_bed.pl < cds_merged_no_centromere.bed)

utr_counts=$(./count_total_length_from_bed.pl < utr_merged.bed)
utr_counts_no_centromere=$(./count_total_length_from_bed.pl < utr_merged_no_centromere.bed)

echo
echo Intron counts = $intron_counts
echo Intron counts without centromeres = $intron_counts_no_centromere
echo
echo Intergenic counts = $intergenic_counts
echo Intergenic counts without centromeres = $intergenic_counts_no_centromere
echo
echo CDS counts = $cds_counts
echo CDS counts without centromeres = $cds_counts_no_centromere
echo
echo UTR counts = $utr_counts
echo UTR counts without centromeres = $utr_counts_no_centromere
echo

# Extract TTAA motif coords across the human genome
#./count_TTAA_motif_from_fasta.pl -f ${GENOME_FASTA} > all_chromosomes_TTAA.bed

# Estimate number of TTA motifs per annotation type excluding centromeric regions
echo "###########################################"
echo Estimating TTAA counts per annotation type
echo "###########################################"

intergenic_ttaa=$(bedtools intersect -a all_chromosomes_TTAA.bed -b intergenic_no_centromere.bed | wc -l)
intronic_ttaa=$(bedtools intersect -a all_chromosomes_TTAA.bed -b introns_merged_no_centromere.bed | wc -l)
utr_ttaa=$(bedtools intersect -a all_chromosomes_TTAA.bed -b utr_merged_no_centromere.bed | wc -l)
cds_ttaa=$(bedtools intersect -a all_chromosomes_TTAA.bed -b cds_merged_no_centromere.bed | wc -l)

echo 
echo Intergenic TTAAs = $intergenic_ttaa
echo Intronic TTAAs = $intronic_ttaa
echo CDS TTAAs = $cds_ttaa
echo UTR TTAAs = $utr_ttaa
echo

echo "#################################"
echo Done!
echo "#################################"




