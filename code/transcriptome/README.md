
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

