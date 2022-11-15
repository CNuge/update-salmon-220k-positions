import gc
import pandas as pd
from seqio import read_fasta, write_fasta
from dna_subset import build_genome_dict
from dna_subset import subset_snp_from_genome, affy_to_major_seq

def build_placed_snp_seqs(snp_data, affy_70mer_data, genome):
    """ iterate across the list of snps
            if placed, then subset a 201mer from the given chromosome
            if unplaced, then pull the 70mer from the affy file
        for each, make nested dict entry of:  {snp : {index:, seq:, major:, minor:, snp_pos:}}
        returns the nested dict, and a list of failed snps
    """
    snp_out_data = {}
    bad_snps = []
    for _, data in snp_data.iterrows():
        if data['CHR'] != 30:
            #subset a 201mer from the given chromosome 
            snp_seq, snp_pos = subset_snp_from_genome( data['CHR'], 
                                                        data['POS'], 
                                                        genome)
            snp_out_data[data['SNP']] = {
                'seq' : snp_seq,
                'index' : data['INDEX'],
                'major' : data['MAJOR'],
                'minor' : data['MINOR'],
                'snp_pos': snp_pos,}
        else:
            #get the 70mer from the affy data
            affy_entry = affy_70mer_data[affy_70mer_data['SNP'] == data['SNP']].iloc[0].to_dict()            
            try:
                snp_seq, major, minor = affy_to_major_seq(affy_entry['Flankk'])
                snp_out_data[data['SNP']] = {
                    'seq' : snp_seq,
                    'index' : data['INDEX'],
                    'major' : major,
                    'minor' : minor,
                    'snp_pos': 36,} #35 leading, snp, 34 trailing is the format
            except:
                bad_snps.append(data)
    return snp_out_data, bad_snps


if __name__ == "__main__":
    #the names, alleles, and locations
    SNP_INPUT = "../data/CIGENE_220K_SNPlocation_majorminor.txt"
    #the full DNA
    GENOME_INPUT = "../data/OLD-GCA_000233375.4_ICSASG_v2_genomic.fna"
    #the 70mers, for unplaced SNPs and double check of the subsetting
    ALT_SNP_DNA = "../data/Axiom_Ssa_220k_Annotation_trimmed.csv"

    snp_data = pd.read_csv(SNP_INPUT, sep = "\t")
    affy_70mer_data = pd.read_csv(ALT_SNP_DNA)

    #7 NAs, just turning them to chr 30 as these have no position info
    snp_data.CHR[snp_data.CHR.isna()] = 30.0
 
    ## 1875 SNPs with no locations, get these from the 70mer file
    # snp_data[snp_data.CHR == 30]

    #will have long headers, and the scaffolds are not appended
    raw_v2_genome_data = read_fasta(GENOME_INPUT)

    v2_genome_cleaned = build_genome_dict(raw_v2_genome_data)
    raw_v2_genome_data = []
    gc.collect()

    snp_out_data_dict, bad_data = build_placed_snp_seqs(snp_data, 
                                                        affy_70mer_data, 
                                                        v2_genome_cleaned)

    #from here, make the header : seq pairs from the snp_out_data_dict
    #and then write to a fasta file for bwa

