import pandas as pd
from seqio import read_fasta


#desc = 'AX-87448985\tSNP_NUMBER 999;MAJOR A;MINOR C;SEQ_POS 101'
def get_poly_pos_from_header(desc):
    """Parse the description in the fasta header and get the alleles and snp position."""
    comps = desc.split(';')

    major = comps[1].split(" ")[-1]
    minor = comps[2].split(" ")[-1]
    pos = comps[3].split(" ")[-1]

    return f"{major}/{minor}", pos

def build_snp_placer_data_frame(snp_seq_info):
    """Parse the fasta info and make the dataframe for snp-placer."""
    out_info = []
    for d in snp_seq_info:
        poly, pos = get_poly_pos_from_header(d['description'])
        outdict = {"SNP" :  d['name'],
                    "Polymorphism" : poly,
                    "bp" : pos,
                    "Sequence" : d['sequence'],}
        out_info.append(outdict)
    return pd.DataFrame(out_info)


#line = ">NC_059442.1 Salmo salar chromosome ssa01, Ssal_v3.1, whole genome shotgun sequence"
def parse_genome_names_v3(NEW_GENOME_NAMES):
    def get_id(line):
        return line.split(" ")[0][1:]

    def get_chr(line):
        if "Salmo salar chromosome" in line:
            chrom = line.split("chromosome ssa")[-1]
            chrom = chrom.split(",")[0]
            return int(chrom)
        #scaffold, assign to chr 30
        return 30

    chr_name_pairs = []
    with open(NEW_GENOME_NAMES) as file:
        for line in file:
            #chromsome
            out = {'ID' : get_id(line), 
                    'chr' : get_chr(line)}
            chr_name_pairs.append(out)

    return pd.DataFrame(chr_name_pairs)


if __name__ == "__main__":

    """
    need the SNP data from the fasta in the format:

    SNP	Polymorphism	bp	Sequence
    CMN1211	G/T	24	TGCATATGGCTCATCACAAATACGCAGAAAAAATGTTGCAGGTGGAGCATCACATGCA
    CMN8988	A/C	51	TGCATATGGCTCTCCTATTCTTTGCCCAGTCATATTCAAGGTTAGAACTCATTTTCTAGGGTTC

    and saved to a .tsv

    also need the header info processed so the NCBI designations can be matched 
    to the 1-29 designations


    """    
    SNP_INPUT = "../data/snp_sequence_records_220k.fasta"
    TSV_OUTPUT = "../data/snp_sequence_records_220k.tsv"
    NEW_GENOME_NAMES_IN = "../data/NEW_genome_chr_names.txt"
    NEW_GENOME_NAMES_OUT = "../data/NEW_genome_chr_ID_number_pairs.txt"

    #read the fasta
    snp_seq_info = read_fasta(SNP_INPUT, return_description = True)
    #build the snp-placer tsv file
    snp_df = build_snp_placer_data_frame(snp_seq_info)
    "write to file"
    snp_df.to_csv(TSV_OUTPUT, sep = "\t", index = False)

    name_pairs = parse_genome_names_v3(NEW_GENOME_NAMES_IN)

    name_pairs.to_csv(NEW_GENOME_NAMES_OUT, sep = "\t", index = False)





