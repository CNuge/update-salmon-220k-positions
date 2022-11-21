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
    merged_map_df = merged_map_df.drop_duplicates()
    #double check for duplicate values onces all together, remove these and move 
    vc_gte_1 = merged_map_df["SNP"].value_counts()
    vc_gte_1 = vc_gte_1[vc_gte_1>1]
    dup_snps = set(vc_gte_1.index)
    map_file_nodup = merged_map_df[[x not in dup_snps for x in merged_map_df["SNP"].values ]]

    map_file_dups = merged_map_df[[x in dup_snps for x in merged_map_df["SNP"].values ]]
    map_file_dups["CHR"] = 30
    map_file_dups["POS"] = 30
    map_file_dups = map_file_dups.drop_duplicates()
    map_file_dups.reset_index(inplace=True, drop = True)
    map_file_dups["POS"] = map_file_dups.index + 220000

    new_map_file_cleaned = pd.concat([map_file_nodup, map_file_dups], ignore_index=True)

    assert len(new_map_file_cleaned) == 220000 #all aboard?

    new_map_file_cleaned.to_csv("../data/Ssa220K_v3_Genome.map", sep = "\t", index = False)

    """ read  in the old map file, do a comparison of the new snp positions to the old
    
        how many moved chromosomes? how many stayed the same? how many unplaced now have homes?"""

    # MAP FORMAT:
    # CHR   SNP ZERO    POS
    #old map in
    old_map = pd.read_csv(OLD_LOCATIONS, sep = "\t")
    old_map.CHR[old_map.CHR.isna()] = 30.0
    old_map["CHR"] = old_map["CHR"].apply(lambda x: int(x))

    #convert new and old maps to dicts and record differences

    new_map_dict = {}
    for i, row in new_map_file_cleaned.iterrows():
        new_map_dict[row["SNP"]] = {"POS": row["POS"], "CHR" : row["CHR"]}

    old_map_dict = {}
    for i, row in old_map.iterrows():
        old_map_dict[row["SNP"]] = {"POS": row["POS"], "CHR" : row["CHR"]}

    total = len(new_map_dict)
    matches = 0
    non_matches = []

    for k, new_v in new_map_dict.items():
        old_v = old_map_dict[k]
        if new_v["CHR"] == old_v["CHR"]:
            matches +=1
        else:
            non_matches.append({"new_CHR" : new_v["CHR"], 
                                 "old_CHR" : old_v["CHR"]})

    matches / total #96.85% go to the same chromosome
    #this seems reasonable