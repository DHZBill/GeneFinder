# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Bill Du

"""
from load import load_seq
dna = load_seq("./data/X73525.fa")

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))



def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide=='A':
    	return 'T'
    if nucleotide=='C':
    	return 'G'
    if nucleotide=='T':
    	return 'A'
    if nucleotide=='G':
    	return 'C'
   


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse=[None]*len(dna)
    for i in range(0,len(dna)):
    	reverse[len(dna)-i-1]=get_complement(dna[i])
    s=''
    reverse=s.join(reverse)
    return reverse

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """    

    for i in range(3,len(dna)-2,3):
    	stop_codon=dna[i:i+3]
    	start_codon=dna[:i]
    	if stop_codon=='TAG' or stop_codon=='TAA' or stop_codon=='TGA':
    		return start_codon
    		break
    	if (i+3) >= len(dna)-2:
        	return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
        The test added below is to make sure that when the given DNA does't start with a start codon,
        and when there's invalid codons in the DNA, the function still gives a correct output.
    >>> find_all_ORFs_oneframe("GGGATGCATGAATGTAGATAGGGGGGGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    

    ORFs=[]
    i=0
    while i<len(dna):
        if dna[i:i+3]=="ATG":
            j=i
            dna2=rest_of_ORF(dna[j:])
           
            if dna2==None:
            	break
            ORFs.append(dna2)
            i=i+len(dna2)
        i=i+3
    return ORFs
    	



def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs_oneframe("GGGATGCATGAATGTAGATAGGGGGGGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORFs=[]
    for i in range (0,3):
    	frame=dna[i:] 
        all_ORFs+=find_all_ORFs_oneframe(frame)
    return all_ORFs
    




    


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    origin=dna
    reverse=str(get_reverse_complement(dna))
    O_ORF=find_all_ORFs(origin)
    R_ORF=find_all_ORFs(reverse)
    Both=O_ORF+R_ORF
    return Both



def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    l = find_all_ORFs_both_strands(dna)
    longest=''
    if len(l)>=1:
	    longest =max(l,key=len)
    return longest



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 
    """
    longest=[]
    for i in range(0,num_trials):
    	shuffled_str=shuffle_string(dna)
    	longest.append(longest_ORF(shuffled_str))
    long_ORF=max(longest,key=len)
    return len(long_ORF)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    protein=''
    for i in range(0,len(dna),3):
	    if dna[i:i+3] in aa_table.keys():
	    	protein += aa_table[dna[i:i+3]]
    return protein


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    l = []
    for i in find_all_ORFs_both_strands(dna):
    	if len(i)>=threshold:
    		l.append(coding_strand_to_AA(i))
    print l
    return l

gene_finder(dna)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
