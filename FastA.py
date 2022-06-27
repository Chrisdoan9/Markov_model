def Read_FastA_sequences(path):
    """This crude function takes as an argument the path to an input file in FastA format that contains
    one or more sequences. The method returns two lists, the first a list of headers (a.k.a. annotation lines), and
    a second corresponding to a list of the sequences themselves.  It is not very smart about spaces, sequence
    position numbers, etc. --- it just nukes all of these
    """
    translation_table = str.maketrans(dict.fromkeys('0123456789* '))
    sequences = []
    headers = []

    header = None
    sequence = ''

    print('Opening', path)

    with open(path, 'r') as f:  # We haven't discussed it in class, but with is a Python context-manager

        line = f.readline()

        while line:

            line = line.strip()

            if line and line[0] == '>':

                if header:
                    headers.append(header)
                    sequences.append(sequence)

                header = line[1:]
                sequence = ''

            else:
                line = line.translate(translation_table)
                sequence += line.upper()

            line = f.readline()

        if header:  # Mop up the last sequence

            headers.append(header)
            sequences.append(sequence)

    print(sequences)
    return headers, sequences


def main():  # Will not execute when FastA is imported as a module

    head, seq = Read_FastA_sequences('seq.pir.txt')

    print(head, seq)


if __name__ == '__main__':
    main()