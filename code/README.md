## README

########################################################################

### Notes on code folder contents: 

## Intro files: 

1) README (this file)
2) packages.R contains all the libraries and packages needed to run scripts in this folder (including the following sub-folders)

## Sub-folders: 

1) genome folder contains all scripts needed for genome-wide analysis
2) transcriptome folder contains all scripts needed for transcriptome analysis
3) comparisons folder contains scripts for regulatory analysis and other combined analyses

#
"genome" folder: 

1) SNPs_filter_maf&cov.R contains the script used to filter the preliminary SNP table for minor allele frequencies and coverage. 
2) indels_filter_maf&cov.R contains the script used to filter the preliminary indel table for minor allele frequencies and coverage.
3) CMHtest_SNPdata.R contains the script used to run the CMH test for the SNP data
4) indels_cmhtest.R contains the script used to run the CMH test for the indel data
5) nuclearSNPs_ManhattanPlot.R is the script used to generate a manhattan plot for the nuclear SNPs
6) mitoSNPs_ManhattanPlot.R is the script used to generate a manhattan plot for the mitochondrial SNPs
7) indels_nuclear_manhattan.R is the script used to generate a manhattan plot for the nuclear indels
8) indels_mito_manhattan.R is the script used to generate a manhattan plot for the mitochondrial indels
9) indels&snps_nuclear_manhattan.R is the script used to generate a manhattan plot that overlays the indels on top of the SNP manhattan plot - incomplete/needs to be fixed.
10) Haplotype_estimation_eNotebook_version.R is the script that was used to estimate the haploype frequencies for each replicate. This script was transferred over after use, and was used iteratively to run for each chromosome and then they were all bound together into a single file. 
11) haplotype_differences.R is the script used to calculate the difference in haplotype frequencies across the genome to use for plotting
12) supp_table_sig_variant_list.R is the script used to create the supplementary table of all significant variants from teh genome data

#
"transcriptome" folder

1) batch_adjust_DESEQ.R is the script that takes the gene count matrix generated in the description below and performs a batch correction based on the sorting date, performs some additional filtering to remove genes with low mapping, and creates a DESEQ object for later analysis
2) DGE_gene_lists.R is the script used to create the significant gene lists for GO-term analysis, and to be used to compare against the genome data
3) heatmap.R is the script used to create the heatmap figures for the transcriptome data
4) rna_supp_table.R is the script used to create the supplementary data of significant transcripts and their distances from implicated genome variants

#
"comparisons" folder: 
1) 00_common_gene_lists.R - this script should be run before running any other scripts in this folder. It takes the significant gene lists from the genomic data and extracts the appropriate annotations for use in later analysis.
2) dna_rna_gene_comparison.R uses files generated from 00_common_gene_lists.R to find genes that are implicated in both genomic and transcriptomic analyses
3) genic_vs_nongenic.R uses files generated from 00_common_gene_lists.R to compare which variants are found within coding vs non-coding regions, and to create a combined list of genic significant variants for genomic and transcriptomic GO-analysis.
4) annotation_comparisons.R uses files generated from 00_common_gene_lists.R to compare the frequency of annotation categories observed in our significant genomic data to the entire annotated genome.
5) annotations_pie_chart.R uses files generated from annotation_comparisons.R to generate a pie chart that visualized the annotations (Figure)
6) FishersTest_annotations.R performs a statistical test on the annotation proportions from annotation_comparisons.R




... to be continued...


#### outputs should be dumped into temp folders

#######################################################################

### Notes on GWAS analysis


########################################################################

### Notes on TWAS Analysis

#### Generation of gene count matrix

#### Notes from CQLS server generation of gene_count_matrix

Used recommendations from Lexogen analysis document, with some mods: 
https://www.lexogen.com/wp-content/uploads/2023/01/015UG009V0271_QuantSeq-3‘-mRNA-Seq_2023-01-24.pdf
(V1 has since been phased out for V2)

##### location of files in shell: /nfs3/IB/Burke_Lab/McHugh/RNASEQ

###### bbduk was used to clean fastq files (*cleaning.sh file in folder on CQLS infrastructure) ###### bbduk is already customized to our server #don't need to unzip files

  SGE_Batch -c "bbduk.sh -Xmx4g in=fastqs/$input out=$output ref=/local/cluster/bbmap/resources/polyA.fa.gz,/local/cluster/bbmap/resources/truseq.fa.gz ftl=12 k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=20 minlength=25" -r bbduk_loop_${prefix}_${i}.sge
  
##### Cleaned files located in filtered_samples folder
Aligned filtered samples to reference genome using STAR (this was the recommended method--also tested HiSat2 and salmon (transcriptome), but STAR had the best alignment)

Used alignment.sh file for alignment in the filtered_samples folder: 

(for sample_type in "OLD" "YOUNG"; do
    for ((i = 1; i <= 12; i++)); do
        sample_name="RNAseq_${sample_type}_rep$(printf "%02d" $i)"
        echo "Processing sample: $sample_name"  # Display the current sample being processed
        SGE_Batch -c "STAR --runThreadN 12 --genomeDir star_genome_ann --outFilterMultimapNmax 20 --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --readFilesIn ${sample_name}_bbduk.fastq.gz --outFileNamePrefix ${sample_name}_ANN_loop. --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMax 100 --alignMatesGapMax 100 --limitBAMsortRAM 1259507435" -r star_${sample_name}.sge
    done
done)

##### output saved to STAR_output folder within filtered_sampes folder
the feature_counts.sh file was used on the alignments, along with the yeast.2023.mRNA.sgd.gff reference file:

for sample_type in "OLD" "YOUNG"; do
    for ((i = 1; i <= 12; i++)); do
        sample_name="${sample_type}_rep$(printf "%02d" $i)"
        echo "Processing sample: $sample_name"  # Display the current sample being processed
        SGE_Batch -c "featureCounts -a yeast.2023.mRNA.sgd.gff -F GFF -t mRNA -g Parent -s 1 -o RNAseq_${sample_name}.ann.featurecounts.stranded.txt RNAseq_${sample_name}_ANN_loop.Aligned.sortedByCoord.out.bam" -r RNAseq_${sample_name}_featcounts.sge
    done
done

rename.sh* used to rename aligned files: for file in *_ANN_loop.Aligned.sortedByCoord.out.bam; do
    mv -- "$file" "${file%_ANN_loop.Aligned.sortedByCoord.out.bam}.bam"
done

##### renamed files were saved within gene_count_matrix directory, featurecounts was used to generate the gene_count_matrix that is used below

featureCounts -F GFF -t mRNA -g Parent -s 1 -a yeast.2023.mRNA.sgd.gff -o gene_count_matrix.txt RNAseq_OLD_rep01.bam RNAseq_OLD_rep02.bam RNAseq_OLD_rep03.bam RNAseq_OLD_rep04.bam RNAseq_OLD_rep05.bam RNAseq_OLD_rep06.bam RNAseq_OLD_rep07.bam RNAseq_OLD_rep08.bam RNAseq_OLD_rep09.bam RNAseq_OLD_rep10.bam RNAseq_OLD_rep11.bam RNAseq_OLD_rep12.bam RNAseq_YOUNG_rep01.bam RNAseq_YOUNG_rep02.bam RNAseq_YOUNG_rep03.bam RNAseq_YOUNG_rep04.bam RNAseq_YOUNG_rep05.bam RNAseq_YOUNG_rep06.bam RNAseq_YOUNG_rep07.bam RNAseq_YOUNG_rep08.bam RNAseq_YOUNG_rep09.bam RNAseq_YOUNG_rep10.bam RNAseq_YOUNG_rep11.bam RNAseq_YOUNG_rep12.bam
echo "  Finished at:           " `date` 

#######################################################################

#### Notes on Regulatory interactions
