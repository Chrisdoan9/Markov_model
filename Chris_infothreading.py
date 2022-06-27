import AA_module
import FastA
import re
import FindInfo


class ThreadProtein(object):
    """ThreadProtein
    VERSION:    Version 5.0 (3 Apr 2020) by Paul Fawcett, but very closely modelled after an earlier PERL program by
                Jeff Elhai.  V4.0 updates to Python 3, V5.0 corrects str conversion of lists in helper methods

    PURPOSE:    Threads the amino acid sequence of one protein # threads?
                through the three-dimensional structure of another

    INPUT FILES:

        Structure file: PDB format. Must have only one chain, corresponding to known sequence.

        Alignment file: FastA format (output as PIR format by Clustal)

            Leading and lagging blanks are OK

            SEQUENCES MUST BE ALIGNED! (i.e. with gaps if necessary)
            First sequence must be that also represented by structure file = Known sequence
            Second sequence must be that to be superimposed on known structure = threaded sequence

            > [header line for sequence with known structure]
            [optional blank line]
            [amino acid sequence, with gaps (-) for alignment]
            * [optional]
            > [header line for sequence to be threaded]
            [optional blank line]
            [amino acid sequence, with gaps (-) for alignment]
            * [optional]

    OUTPUT FILES: Description of format of output files

        New structure file: PDB format

            NO CHANGE: Amino acids common between the two sequences go in chain A
            REPLACEMENTS: Amino acids changed in threaded sequence go in chain B
                Amino acid called new name in first atom of residue (others unchanged)
            INSERTIONS: Amino acids in threaded sequence but not original go in chain C
                No attempt made to find positions of the new atoms.
                Only first atom given
                Inserted atom given meaningless number
            DELETIONS: Amino acids in original but not in threaded sequence go in chain D

            Informational lines (e.g. Title, source, etc) copied unchanged
            therefore information may no longer be appropriate!!
    """

    translation_table = str.maketrans(dict.fromkeys('-', ''))

    def __init__(self, structure_file=r'1ei1.pdb',
                 alignment_file=r'post_align.txt',
                 output_file=r'superimpose_DNA_Gyrase.pdb'):

        self.structure_file = structure_file
        self.alignment_file = alignment_file
        self.output_file = output_file

        headers, sequences = FastA.Read_FastA_sequences(alignment_file)

        self.known_seq_header = headers[0]
        self.threaded_seq_header = headers[1]

        self.known_seq = sequences[0]
        self.threaded_seq = sequences[1]

        # mutan_position = set[(198, 10)]
        # for pos in mutan_position:
        #     self.__translate_positions(pos)
        #     print(self.comparison_summary)

        self.comparison_summary = self.__analyze_alignment()

        # self.mutation_indeces = self.__read_mutations()
        self.__Modify_comparison_summary()
        self.__create_new_structure_file()
        print('Finished writing new PDB file', self.output_file)

    def __analyze_alignment(self):

        """ Analyze alignment at each position
        Put I (insertion) where gap in sequence with known structure
        Put D (delection) where gap in sequence to be threaded
        Put M (match) where both sequences match
        Put R (replacement) where sequences differ

        relies on instance variable self.known_seq and self.threaded_seq
        """

        comparison_summary = ''
        gap_char = '-'

        print('\nAnalyzing alignment:\n')
        print(self.known_seq)
        print(self.threaded_seq)

        if not len(self.known_seq) == len(self.threaded_seq):
            quit("The two sequences must be the same length. Don't forget to align them!")

        # print('Length of input sequences is', len(self.known_seq))

        # First we create a string that summarizes the alignment in a position-by-position fashion

        for i in range(len(self.known_seq)):

            known_char = self.known_seq[i]
            threaded_char = self.threaded_seq[i]

            if known_char == gap_char:
                comparison_summary += 'I'
            elif threaded_char == gap_char:
                comparison_summary += 'D'
            elif known_char == threaded_char:
                comparison_summary += 'M'
            else:
                comparison_summary += 'R'

        comparison_summary += 'M'  # Adds a final match so as to correctly handle the termination sequence
        # self.__read_mutations() # Q added

        self.comparison_summary = comparison_summary
        print('Comparison summary before')
        print(self.comparison_summary)
        # print(self.__translate_positions())

        # quit()
        return comparison_summary

    def __Modify_comparison_summary(self):
        _, _, _, indices = FindInfo.Calculate_info_and_consensus(input_filepath='Chris_Gyrb_aligned.fa.txt',
                                                                 alphabet='ACDEFGHIKLMNPQRSTVWY')
        count_del = 0  # to find equivalent position of conserved residues in new alignment file.
        while self.threaded_seq[count_del] == '-':
            count_del = count_del + 1
        print('del number', count_del)
        for i in range(len(indices)):
            indices[i] = indices[i] + count_del

        print('Indices conserved amino acids', indices)
        self.translated = indices

        # for pos in self.mutation_indeces:
        #     self.translated.append(self.__translate_positions(pos))
        # print('here', self.__translate_positions(1))
        #
        temp = list(self.comparison_summary)
        for i in self.translated:
            temp[i] = 'X'

        self.comparison_summary = ''.join(temp)
        print('Comparison summary after')
        print(self.comparison_summary)
        return self.comparison_summary

    def __create_new_structure_file(self):  # change 4 letters on header before upload to jmol

        # gap_char = '-'
        gaps_reached_in_known_seq = 0
        terminus_found = False
        chains_printed = False
        previous_line = ''

        last_residue = len(self.known_seq.translate(self.translation_table))
        atom_lines = dict()

        atom_lines['A'] = []
        atom_lines['B'] = []
        atom_lines['C'] = []
        atom_lines['D'] = []
        atom_lines['E'] = []

        with open(self.structure_file, 'r') as old_structure:
            with open(self.output_file, 'w') as new_structure:

                line = old_structure.readline()

                while line:

                    line = line.strip()
                    line_type = line[0:4]

                    if line_type == 'ATOM' or line_type == 'TER ':
                        residue = int(line[22:26]) - 1 + gaps_reached_in_known_seq  # Python residues begin at 0

                        if line_type == 'TER ':  # PDB file numbers terminal residue
                            residue += 1  # same as last residue. Fix this.
                            terminus_found = True

                        # print(residue, 'is', self.comparison_summary[residue])
                        residue_status = self.comparison_summary[residue]

                        if residue_status == 'M':

                            print(line, '\n', file=new_structure)
                            previous_line = line
                            line = old_structure.readline()

                        elif residue_status == 'R':

                            new_residue = self.threaded_seq[residue]

                            line = self.__update_residue_name(line, new_residue)
                            line = self.__update_chain_id(line, 'B')
                            atom_lines['B'].append(line)
                            previous_line = line
                            line = old_structure.readline()

                        elif residue_status == 'X':

                            new_residue = self.threaded_seq[residue]

                            line = self.__update_residue_name(line, new_residue)
                            line = self.__update_chain_id(line, 'E')
                            atom_lines['E'].append(line)
                            previous_line = line
                            line = old_structure.readline()

                        elif residue_status == 'D':
                            newline = self.__update_chain_id(line, 'D')
                            atom_lines['D'].append(newline)
                            previous_line = line
                            line = old_structure.readline()

                        elif residue_status == 'I':
                            last_residue += 1
                            new_residue = self.threaded_seq[residue]
                            if line_type == 'TER ':
                                newline = previous_line
                            else:
                                newline = line

                            newline = self.__update_residue_name(newline, new_residue)
                            newline = self.__update_chain_id(newline, 'C')
                            newline = self.__update_residue_number(newline, last_residue)
                            gaps_reached_in_known_seq += 1
                            atom_lines['C'].append(newline)

                        else:
                            print('Error in:\n{}\n\tResidue {}\n\tResidue status: {}'.format(line, residue,
                                                                                             residue_status))

                            line = old_structure.readline()

                    else:
                        if terminus_found and not chains_printed:
                            self.__output_extra_chains(atom_lines, new_structure)
                            chains_printed = True
                        print(line, file=new_structure)
                        line = old_structure.readline()

    @staticmethod
    def __update_residue_name(line, residue_name):

        line = list(line)  # string type doesn't support item assignment, so convert to list

        if residue_name == 'INS' or residue_name == 'DEL':
            line[17:20] = residue_name
        else:
            line[17:20] = AA_module.One_letter_to_three_letter_code(residue_name)

        return ''.join(line)

    @staticmethod
    def __update_chain_id(line, chain_id):

        line = list(line)  # string type doesn't support item assignment, so convert to list

        if not len(chain_id) == 1:  # Chain ID should ALWAYS be length 1
            print('ERROR: Chain ID {} not of length 1'.format(chain_id))
            quit()
        line[21] = chain_id

        return ''.join(line)

    @staticmethod
    def __update_residue_number(line, residue_number):

        line = list(line)  # string type doesn't support item assignment, so convert to list

        new_number = "%4d" % residue_number
        line[22:26] = new_number
        return ''.join(line)

    @staticmethod
    def __output_extra_chains(atom_lines, filehandle):

        # print('Outputing remaining chains')

        for chain in sorted(atom_lines.keys()):

            for atom in range(len(atom_lines[chain]) - 1):
                print(atom_lines[chain][atom], file=filehandle)

    def __translate_positions(self, pos):
        com_sum_pos = 0
        mloti_pos = 0
        while mloti_pos < int(pos):
            if self.comparison_summary[com_sum_pos] != 'D':
                mloti_pos += 1
            # print(com_sum_pos, mloti_pos)
            com_sum_pos += 1
        return com_sum_pos - 1  # -1

        # dict = {}
        # #number_D = 0
        # d = 1
        # for i in range(len(self.comparison_summary)):
        #     if i != 'D':
        #         dict[i] = dict.values[i-d]
        #     if i == 'D':
        #         #number_D = number_D + 1
        #         d = d - 1

        # return dict

    # def translate(self, cs, ml_pos):
    #
    #
    #     pass
    # something here
    def __read_mutations(self):
        input_file = open('UDPGD-mutants.txt', 'r')
        regex = re.compile(r'\w{3}\s+(\d+)\s\w{3}')
        to_Search = input_file.read()
        # print(set(regex.findall((to_Search))))
        return set(regex.findall((to_Search)))


def main():
    ThreadProtein()  # Just use the default filenames for now


if __name__ == '__main__':
    main()

# write program convert the number from regrex to position
# Paul uploaded edited file but it didn't read the mutant file!
