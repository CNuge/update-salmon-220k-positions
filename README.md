# update-salmon-220k-positions
## Orient the SNPs from the Atlantic salmon 220K chip in v3 of the Atlantic salmon genome

## The plan:

[x] - Set up repository on compute canada

[x] - get copy of the v3 Atlantic salmon genome downloaded and added to /data/ folder
        -`NEW-GCF_905237065.1_Ssal_v3.1_genomic.fna`

[x] - add copy of the v2 (current) genome to the /data/ folder
        -`OLD-GCA_000233375.4_ICSASG_v2_genomic.fna`

[x] - add copy of the .map file for the current chip positions to the /data/ folder
    - need a completely unfiltered map file, with the full list of 220k SNPs.
    - `CIGENE_220K_SNPlocation_majorminor.txt` Even better, has the relevant info for the alleles so no ped file needed!
    - COLUMNS: POS MAJOR   MINOR   CHR INDEX   SNP

[] - build a key / dictonary file with the names of the chromosomes from the three files, make sure they're matched unambigiously
    OLD_genome_chr_names.txt, NEW_genome_chr_names.txt, and unique of `CIGENE_220K_SNPlocation_majorminor.txt` CHR


[x] - write a program to create a new fasta file from the old genome file and the map file
        - this step will need some programming 
        - take 100bp up and downstream of the exact SNP position
        - make a header with : >"SNP_NAME";"MAJOR";"MINOR"
        - EDGE case - old scaffolds, append them together linearly to make a singe CHR30 and work off the POS information from there.
            - not positioned, therefore I'm going off the 70mers from the affy file for these SNPs.
        - Will need to handle the new scaffolds as well, placing things onto CHR 30
            -probably just want a function to do this.

[] - spot checks to make sure the SNPs are where the fasta says they are
    - make sure major allele is sitting at the 101st or 36th base pair depending on the two categories

[] - run bwa (on sharcnet), aligning the fasta file from previous step to the v3 genome

[] - use SNP-placer to get the updated position file in .vcf format.

[] - take the relevant information from the SNP-placer outputs vcfs and make a new map file
    - as part of this, make a new 'exclude_v3_genome.tsv' file that can be used to pull SNPs from the PED files using plink 

[] - generate some summary stats regarding the number of ambigious and unambigous SNP positions in the new genome
        - how many SNPs are lost, of these how many were due to no placement and how many were due to multiple placements

[] - from here, carry the new information through to the salmon pop gen work, redoing the manhattan plots for the Fst, the LD calculations, and anything else required.
    -need to double check the orders of the SNPs, sort the final map file by their original positions (so that no weird bugs get introduced)


## Things to test for / keep an eye on

- confirm that the fasta files have the SNPs in the 101st (of 201) bp positions. Spot check a few and make sure the major alleles are there.
    - check against the file: `Axiom_Ssa_220k_Annotation_trimmed.csv` which appears to have the affymetrix primers for the SNPs

- when you have the final bp positions in the new genome, spot check the file to make sure that the alleles for the SNP are in fact found at the specific bp indicated. 

- make sure the chromosome names are matched on the front and back end correctly, keep an eye out for any weird changes or things that are out of order in the files.
    -  don't assume positional equivalency, get eyes on all the name matches

- *** IF THE ALIGNMENT IS A REVERSE COMPLIMENT, MAKE SURE YOU SWITCH THE MAJOR/MINOR ALLELES TO THEIR REVERSE COMPLIMENTS WHEN YOU GENERATE THE OUTPUT
    - want the sequences to stay consistent.