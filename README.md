# update-salmon-220k-positions
## Orient the SNPs from the Atlantic salmon 220K chip in v3 of the Atlantic salmon genome

## The plan

[] - Set up repository on compute canada

[] - get copy of the v3 Atlantic salmon genome downloaded and added to /data/ folder

[] - add copy of the v2 (current) genome to the /data/ folder

[] - add copy of the .map file for the current chip positions to the /data/ folder

[] - write a program to create a new fasta file from the old genome file and the map file
        - this step will need some programming 
        - take 100bp up and downstream of the exact SNP position
        - make a header with : >"SNP_NAME";"MAJOR";"MINOR"

[] - run bwa (on sharcnet), aligning the fasta file from previous step to the v3 genome

[] - use SNP-placer to get the updated position file in .vcf format.

[] - take the relevant information from the SNP-placer outputs vcfs and make a new map file
    - as part of this, make a new 'exclude_v3_genome.tsv' file that can be used to pull SNPs from the PED files using plink 

[] - generate some summary stats regarding the number of ambigious and unambigous SNP positions in the new genome
        - how many SNPs are lost, of these how many were due to no placement and how many were due to multiple placements

[] - from here, carry the new information through to the salmon pop gen work, redoing the manhattan plots for the Fst, the LD calculations, and anything else required.
    -need to double check the orders of the SNPs, sort the final map file by their original positions (so that no weird bugs get introduced)