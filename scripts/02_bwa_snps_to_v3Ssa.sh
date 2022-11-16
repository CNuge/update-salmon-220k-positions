#!/bin/bash
#SBATCH --time=0-03:00 # time (DD-HH:MM)
#SBATCH --job-name=bwa_mem_v_genome
#SBATCH --output=logs/bwa_mem_v_genome_%J.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=125G

#add the scripts folder to the search path
export PATH="$PATH:$(pwd)/bin"

echo "loading necessary modules"
module load StdEnv/2020  
module load gcc/9.3.0   
module load bwa/0.7.17 
module load samtools/1.13


#location of the reference genome to align to
#your genome file here!
genome_file="../data/NEW-GCF_905237065.1_Ssal_v3.1_genomic.fna"
#input SNP data
snp_file="../data/snp_sequence_records_220k.fasta"
cores=32 #make this equal to the cpus-per-task argument given to sbatch on line 5
#location where the alignment outputs should be saved to (in sorted bam format)
#by default to a subfolder in original data directory with name bamfiles/
out_folder="../data/alignment_outputs/"
out_file="220ksnps_align_v3_ssa_genome"

echo "path to the genome file: $genome_file";
echo "path to the snp fasta file to be processed: $snp_file";
echo "number of cores being used: $cores";
echo "writing outputs to folder: $out_folder"

#build the output folder if non existant
if [ ! -d $outfolder ]; then
  mkdir -p $outfolder;
fi

#index the reference file
#bwa index $genome_file

sam_outfile=$outfolder$out_file".sam";	 # the name of the .bam output
sam_sorted_outfile=$outfolder$out_file".sorted.sam";	 # the name of the .bam output

echo "saving data to file:" $sam_outfile;

#burrows wheeler mem alignment to reference genome and output sorting with samtools
bwa mem -t $cores $genome_file $snp_file > $sam_outfile

samtools sort $sam_outfile -o $sam_sorted_outfile -T $f -@ $cores -m 3G;

#bwa mem:
#-t #threads
#samtools:
#@T gives number of threads, -m 3G gives the maximum mempry per thread (currently hardcoded, increase)

