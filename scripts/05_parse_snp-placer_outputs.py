import pandas as pd




if __name__ == "__main__":

    #info to populate the new map file
    NEW_VCF_POS = "../data/v3_genome_locations.vcf"

    #use this for the comparison of locations
    #also to grab the dropped SNPs (add to 30)
    SNP_SEQ_RECORDS = "../data/snp_sequence_records_220k.tsv"

    #use this to update the map file to #s
    NEW_CHR_NAMES = "../data/NEW_genome_chr_ID_number_pairs.txt"

    #make this file
    NEW_MAP_OUT = "../data/v3_genome_snps.map"

    """ read in the vcf, convert to a map file"""

    v3_vcf_df = pd.read_csv(NEW_VCF_POS, sep = '\t')

    """ read in the list of SNPs (input to snp-placer), find those that were not on the chromosomes and 
        add to the end of the new map file as chr 30"""


    """ read  in the old map file, do a comparison of the new snp positions to the old
    
        how many moved chromosomes? how many stayed the same? how many unplaced now have homes?"""