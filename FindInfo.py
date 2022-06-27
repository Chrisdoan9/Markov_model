from FastA_V2 import FastA
from math import log


def Calculate_info_and_consensus(input_filepath='Chris_Gyrb_aligned.fa.txt', alphabet='ACDEFGHIKLMNPQRSTVWY', threshold = 4):

    gap = '-'                           # Gap characters will need to be treated specially
    consensus = []                      # The consensus sequence of the alignment will eventually be here
    gapped_consensus = []               # A gapped version of the alignment's consensus
    sequence_matrix = []                # A list of the sequences we'll be analyzing
    information = []                    # A position-by-position information value deduced from the alignment

    sequences = FastA(input_filepath)   # Instantiate a FastA object that we can iterate across to grab sequences

    for annotation, sequence in sequences:      # Now actually grab the sequences

        sequence_matrix.append(sequence)        # Stash them all in memory


    alignment_length = len(sequence_matrix[0])
    # We're just assuming here that all the sequences are properly aligned and of the same length. No error-checking!

    for position in range(alignment_length):           # Analyze the alignment on a position-by-position basis

        symbol_counts = {symbol: 0 for symbol in alphabet}    # A dict comprehension to initialize all counts to zero
        symbol_sums = 0                                       # A count of how often each non-gap symbol is encountered
        gap_count = 0                                         # A count of how often a gap symbol is encountered

        for sequence in sequence_matrix:                # At each position, we must consider each sequence in turn
                                                        # i.e. we are analyzing a column of the alignment
            current_symbol = sequence[position]         # what symbol are we looking at right now?

            if current_symbol != gap:

                symbol_counts[current_symbol] += 1      # Count both symbol-specific and total number of symbols
                symbol_sums += 1

            else:
                gap_count += 1

        information.append(log(len(alphabet), 2))

        # We will start with the maximum uncertainty, and subtract uncertainty resulting from our observations as we go

        # Note that in the following we are essentially treating gap characters as pseudocounts, as they should
        # logically be contributing to the uncertainty.  The "gap pseudocounts" are divvied up equally between all
        # of the symbols in our alphabet.

        for symbol in alphabet:

            symbol_freq = (symbol_counts[symbol] + (gap_count / len(alphabet))) / (symbol_sums + gap_count)



            try:

                information[position] += symbol_freq * log(symbol_freq, 2)



            except ValueError:
            # The ValueError occurs in the case when the symbol_freq is zero -- taking the log of zero
            # is a math domain error.  This will only happen both when a particular symbol was never observed in
            # that column of the alignment AND there were no gaps to contribute a gap pseudocount.
                pass        # ..so we can just disregard this iteration, as it doesn't effect the overall information

        consensus.append(max(symbol_counts, key=lambda k: symbol_counts[k]))


        # This clever little lambda function will return the key associated with the maximum value found in a dict
        # This will give us the symbol associated with the highest count at the current position - the consensus!
        # It's a trick worth remembering. Lambda functions were beyond the scope of BNFO601, but can be quite handy.
        # They are like little nameless one-line functions.

        if gap_count > max(symbol_counts.values()):
            gapped_consensus.append(gap)
        else:
            gapped_consensus.append(max(symbol_counts, key=lambda k: symbol_counts[k]))

    consensus_string = ''.join(consensus)       # stringify
    gapped_consensus_string = ''.join(gapped_consensus)
    print('\n', gapped_consensus_string)
    print(consensus_string)
    print(information)
    indice_conserved=[]
    # threshold = 4
    for i in range(len(information)):
        if information[i] > threshold:
            indice_conserved.append(i)
    print(indice_conserved)
    # print(consensus_string[74])
    return consensus_string, gapped_consensus_string, information, indice_conserved


if __name__ == '__main__':

    Calculate_info_and_consensus()