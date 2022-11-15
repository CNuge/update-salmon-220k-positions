import gc
import pandas as pd
from seqio import read_fasta, write_fasta


#header = ">CM003279.1 Salmo salar isolate Sally breed double haploid chromosome ssa01, whole genome shotgun sequence"
def subset_v2_chr_header(header):
    """Take a long form header from the chr 2 genome 
        and subset just the ssa style chr designation."""
    front = header.split(',')[0]
    chrom = front.split("haploid chromosome ")[-1]

    if chrom[:3] != "ssa":
        raise ValueError(f"Warning header subset does not conform to chr"+\
                            "naming convention: {chrom}")

    return chr


def build_subset_dict(genome_dict, merged_scaffold_name= "ssa30"):
    """Take a dictonary of fasta entries, simplify the chromosome names, and 
        merge all of the unplaced scaffolds (in order) to make an additional pseudochromosome"""
    subset_dict = {}
    for k, v in genome_dict.items():
        if "chromosome ssa" in k:
            new_k = subset_v2_chr_header(k)
            subset_dict[new_k] = v

        elif merged_scaffold_name in subset_dict.keys():
            subset_dict[merged_scaffold_name] += v
            
        else:
            subset_dict[merged_scaffold_name] = v

    return subset_dict


if __name__ == "__main__":

    SNP_INPUT = "../data/CIGENE_220K_SNPlocation_majorminor.txt"
    GENOME_INPUT = "../data/OLD-GCA_000233375.4_ICSASG_v2_genomic.fna"

    snp_data = pd.read_csv(SNP_INPUT, sep = "\t")

    #will have long headers, and the scaffolds are not appended
    raw_v2_genome_data = read_fasta(GENOME_INPUT)

    v2_genome_cleaned = build_subset_dict(raw_v2_genome_data)
    raw_v2_genome_data = []
    gc.collect()

    