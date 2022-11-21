import pandas as pd
from seqio import read_fasta



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
    
    #headers can be used to get the info for the output df
    snp_seq_info = read_fasta(SNP_INPUT, return_description = True)

    



