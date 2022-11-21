import pandas as pd




if __name__ == "__main__":

    #info to populate the new map file
    NEW_VCF_POS = "../data/v3_genome_locations.vcf"

    #use this to grab the dropped SNPs (add to 30)
    SNP_SEQ_RECORDS = "../data/snp_sequence_records_220k.tsv"

    #use this for the comparison of locations
    OLD_LOCATIONS = "../data/CIGENE_220K_SNPlocation_majorminor.txt"

    #use this to update the map file to #s
    NEW_CHR_NAMES = "../data/NEW_genome_chr_ID_number_pairs.txt"

    #make this file
    NEW_MAP_OUT = "../data/v3_genome_snps.map"

    #read in the chr name map
    v3_chr_name_df = pd.read_csv(NEW_CHR_NAMES, sep = '\t')
    v3_chr_name_map = dict(zip(v3_chr_name_df["ID"].values, v3_chr_name_df["chr"].values))

    """ read in the vcf, convert to a map file"""
    v3_vcf_df = pd.read_csv(NEW_VCF_POS, sep = "\t")
    v3_vcf_df.columns
    #pull the IDS
    v3_vcf_df["SNP"] = v3_vcf_df["ID"].apply(lambda x : x.split("_")[0])
    v3_vcf_df["CHR"] = v3_vcf_df["#CHROM"].apply(lambda x: v3_chr_name_map[x])
    v3_vcf_df["ZERO"] = 0

    map_file_output = v3_vcf_df[["CHR", "SNP", "ZERO", "POS"]]

    """ read in the list of SNPs (input to snp-placer), find those that were not on the chromosomes and 
        add to the end of the new map file as chr 30"""
    snp_seq_df = pd.read_csv(SNP_SEQ_RECORDS, sep = "\t")

    mapped_snps = set(v3_vcf_df['SNP'].values)
    #123 unmapped SNPS
    unmapped_snps = snp_seq_df[[x not in mapped_snps for x in snp_seq_df['SNP'].values]]
    unmapped_snps.reset_index(inplace=True)
    unmapped_snps['CHR'] = 30
    unmapped_snps['POS'] = unmapped_snps.index
    unmapped_snps["ZERO"] = 0
    map_file_addendum = unmapped_snps[["CHR", "SNP", "ZERO", "POS"]]

    merged_map_df = pd.concat([map_file_output, map_file_addendum])
    #double check for duplicate values onces all together, remove these and move 
    vc_gte_1 = merged_map_df["SNP"].value_counts()
    vc_gte_1 = vc_gte_1[vc_gte_1>1]
    vc_gte_1



    """ read  in the old map file, do a comparison of the new snp positions to the old
    
        how many moved chromosomes? how many stayed the same? how many unplaced now have homes?"""

    # MAP FORMAT:
    # CHR   SNP ZERO    POS