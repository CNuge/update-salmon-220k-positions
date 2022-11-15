import re

#header = "Salmo salar isolate Sally breed double haploid chromosome ssa01, whole genome shotgun sequence"
def subset_v2_chr_header(header):
    """Take a long form header from the chr 2 genome 
        and subset just the ssa style chr designation."""
    front = header.split(',')[0]
    chrom = front.split("haploid chromosome ")[-1]

    if chrom[:3] != "ssa":
        raise ValueError(f"Warning header subset does not conform to chr"+\
                            "naming convention: {chrom}")

    return chrom

def build_genome_dict(genome_list):
    """Take a dictonary of fasta entries, simplify the chromosome names.
        Drops the unnplaced scaffolds"""
    subset_dict = {}

    for x in genome_list:
        desc = x['description']
        if "chromosome ssa" in desc:
            seq = x['sequence']
            new_k = subset_v2_chr_header(desc)
            subset_dict[new_k] = seq
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


"""
## Tests for the slicing in the subset genome functions

#normal case
ex_chr = "A" * 200 + "T" * 100 + "C" + "T" * 100 + "A" * 200
pos = 300 #301 - 1
ex_chr[pos]

front_edge = pos - 100
back_edge = pos + 101

assert len(ex_chr[front_edge:back_edge]) == 201
assert ex_chr[front_edge:back_edge][100] == "C"

#short front
front_ex_chr = "T" * 60 + "C" + "T" * 100 + "A" * 200
pos = 60 #61-1
front_edge = 0 
back_edge = pos + 101

front_ex_chr[front_edge:back_edge][pos]
assert front_ex_chr[front_edge:back_edge][pos] == "C"
assert len(front_ex_chr[front_edge:back_edge]) == 161

#short back
back_ex_chr = "A" * 200 + "T" * 100 +  "C" + "T" * 55
pos = 300

front_edge = pos - 100
back_edge = len(back_ex_chr)

assert back_ex_chr[front_edge:back_edge][100] == "C"
assert len(back_ex_chr[front_edge:back_edge]) == 156
"""

def subset_snp_from_genome(chr, pos, genome, flank = 100):
    """Take a SNP position and a given chromosome 
        and subset a the snp with flanks of the given size.
        
        Input position is assumed to be 1 indexed (and is converted to 0 index)"""
    chr_key = f"ssa{int(chr)}"
    chr_seq = genome[chr_key]
    adj_pos = pos - 1
    #literal edge cases - early in chr
    if adj_pos < flank :
        front_edge = 0
        back_edge = adj_pos + flank + 1
        snp_pos = adj_pos + 1 #converts it back to 1 index for output
    #literal edge cases - late in chr
    elif adj_pos >= (len(chr_seq) - flank):
        front_edge = pos - flank
        back_edge = len(chr_seq)
        snp_pos = flank + 1
    #normal case - plenty of flank
    else:
        front_edge = adj_pos - flank
        back_edge = adj_pos + flank + 1
        snp_pos = flank + 1
    snp_seq = chr_seq[front_edge:back_edge]
    return snp_seq, snp_pos
