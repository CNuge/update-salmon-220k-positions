#! /bin/bash

#command line inputs to run snp-placer



cd /mnt/c/Users/camnu/bin/snp-placer

python filter_sam_file.py ../update-salmon-220k-positions/data/alignment_outputs/220ksnps_align_v3_ssa_genome.sorted.sam

mv one_location_alignments.sam ../update-salmon-220k-positions/data/alignment_outputs
mv unaccounted_alignments.sam ../update-salmon-220k-positions/data/alignment_outputs

cd /mnt/c/Users/camnu/bin/update-salmon-220k-positions/data/alignment_outputs
#drop duplicates, seems 2x of each entry in file, not sure why
cat -n one_location_alignments.sam | sort -uk2 | sort -nk1 | cut -f2- > one_location_alignments.dedup.sam
cat -n unaccounted_alignments.sam | sort -uk2 | sort -nk1 | cut -f2- > unaccounted_alignments.dedup.sam

cd /mnt/c/Users/camnu/bin/snp-placer

python place_snps.py -s ../update-salmon-220k-positions/data/alignment_outputs/one_location_alignments.sam -p ../update-salmon-220k-positions/data/snp_sequence_records_220k.tsv

mv placed_snps.vcf ../update-salmon-220k-positions/data/v3_genome_locations.vcf