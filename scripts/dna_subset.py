

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

def build_genome_dict(genome_dict):
    """Take a dictonary of fasta entries, simplify the chromosome names.
        Drops the unnplaced scaffolds"""
    subset_dict = {}
    for k, v in genome_dict.items():
        if "chromosome ssa" in k:
            new_k = subset_v2_chr_header(k)
            subset_dict[new_k] = v
    return subset_dict

# omitted as there are no exactly placed snps from the scaffolds
# could use this to build a pseudochromosome of the scaffolds later
#        elif merged_scaffold_name in subset_dict.keys():
#            subset_dict[merged_scaffold_name] += v           
#        else:
#            subset_dict[merged_scaffold_name] = v


#seq = 'AGATCAAGGGTCCAGTGAGAGATCAGACGTGTGAC[C/T]GGAAACTGGAAACTTACACTTCTGAGGAGGGGGG'
def affy_to_major_seq(seq):
    """Take in an affymetrix formatted 70mer and convert to a standard DNA string format.
        Returns the new string and the associated major and minor alleles.
        
        First allele is minor, second is major."""

    #split on: [ ] and /
    s1, minor, major, s2 = res = re.split('\[|\]|\/', seq)

    outseq = s1+major+s2

    return outseq, major, minor
