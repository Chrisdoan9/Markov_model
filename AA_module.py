def One_letter_to_DNA_codons(amino_acid):
    """ Take as an argument a symbol corresponding to a one letter amino acid code and return a tuple
     corresponding to all of the DNA codons corresponding to that amino acid in the standard genetic code
    """

    codons = {
        'I': ('ATT', 'ATC', 'ATA'),
        'L': ('CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'),
        'V': ('GTT', 'GTC', 'GTA', 'GTG'),
        'F': ('TTT', 'TTC'),
        'M': ('ATG',),
        'C': ('TGT', 'TGC'),
        'A': ('GCT', 'GCC', 'GCA', 'GCG'),
        'G': ('GGT', 'GGC', 'GGA', 'GGG'),
        'P': ('CCT', 'CCC', 'CCA', 'CCG'),
        'T': ('ACT', 'ACC', 'ACA', 'ACG'),
        'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
        'Y': ('TAT', 'TAC'),
        'W': ('TGG',),
        'Q': ('CAA', 'CAG'),
        'N': ('AAT', 'AAC'),
        'H': ('CAT', 'CAC'),
        'E': ('GAA', 'GAG'),
        'D': ('GAT', 'GAC'),
        'K': ('AAA', 'AAG'),
        'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
        '*': ('TAA', 'TAG', 'TGA')}

    return codons[amino_acid.upper()]


def One_letter_to_three_letter_code(amino_acid):
    """ Take as an argument a symbol corresponding to a one letter amino acid code and return a string
     corresponding that amino acid in the three-letter system of amino acid nomenclature.
    """

    one_letter = \
        {'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN',
            'D': 'ASP', 'N': 'ASN', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR',
            'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET', 'A': 'ALA',
            'G': 'GLY', 'P': 'PRO', 'C': 'CYS', '*': '*'}

    return one_letter[amino_acid.upper()]


def Three_letter_to_one_letter_code(amino_acid):
    """ Take as an argument a string corresponding that amino acid in the three-letter system of amino acid
    nomenclature and return a string corresponding to a one letter amino acid code
    """

    three_letter = \
        {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', '*': '*'}

    return three_letter[amino_acid.upper()]

def Report_molecular_weight_of_protein(protein_sequence):
    """
    >>> print(Report_molecular_weight_of_protein('ARNDCEQGHILKMFPSTWYV'))
    2738.0999999999995
    """
    pass  # You'll need to provide the stuffing here

if __name__ == '__main__':

    import doctest
    doctest.testmod()
    print(Report_molecular_weight_of_protein('ARNDCEQGHILKMFPSTWYV'))

# Molecular Weights of the Amino Acids
# Amino Acid	3-letter Code	1-letter Code	Molecular Weight (g/mol)
# Alanine	Ala	A	89.1
# Arginine	Arg	R	174.2
# Asparagine	Asn	N	132.1
# Aspartate	Asp	D	133.1
# Cysteine	Cys	C	121.2
# Glutamate	Glu	E	147.1
# Glutamine	Gln	Q	146.2
# Glycine	Gly	G	75.1
# Histidine	His	H	155.2
# Isoleucine	Ile	I	131.2
# Leucine	Leu	L	131.2
# Lysine	Lys	K	146.2
# Methionine	Met	M	149.2
# Phenylalanine	Phe	F	165.2
# Proline	Pro	P	115.1
# Serine	Ser	S	105.1
# Threonine	Thr	T	119.1
# Tryptophan	Trp	W	204.2
# Tyrosine	Tyr	Y	181.2
# Valine	Val	V	117.1