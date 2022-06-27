__author__ = 'Wombat'

import math
import TMM_HMM_dicts
import FastA_V2
import re
from pickle import load


class HMM(object):

    def __init__(self, sequence=None, states=None, emissions=None):
        print("Loading model file", 'acid_dict')  # De-pickle our model dict, comment out to run with TMM_HMM_dicts.py
        with open('acid_dict', 'rb') as f:
            model = load(f)
        if sequence:
            self.observed_sequence = sequence.upper()  # Could optionally do some cleaning-up here of the sequence
        else:
            self.observed_sequence = []

        self.observed_sequence_length = len(self.observed_sequence)

        if states:
            # self.states = TMM_HMM_dicts.states
            self.states = model
            self.log_states = {outer_k: {inner_k: math.log(inner_v) for (inner_k, inner_v) in outer_v.items()} for
                               (outer_k, outer_v) in self.states.items()}  # nested dict comprehension to log-transform

        else:
            # The default here is a die-rolling occasionally dishonest casino
            self.states = {
                "S": {
                    "F": 0.5,
                    "L": 0.5,
                },
                "F": {
                    "F": 0.95,
                    "L": 0.05,
                },
                "L": {
                    "L": 0.90,
                    "F": 0.10,
                }
            }

        if emissions:
            self.emissions = TMM_HMM_dicts.emissions
            self.log_emissions = {outer_k: {inner_k: math.log(inner_v) for (inner_k, inner_v) in outer_v.items()} for
                                  (outer_k, outer_v) in self.emissions.items()}

        else:
            #  Again, default to a die-rolling occasionally dishonest casino
            self.emissions = {

                "S":  # 'F' indicates a fair die
                    {
                        "_": 1
                    },
                "F":
                    {
                        "1": 1 / 6,
                        "2": 1 / 6,
                        "3": 1 / 6,
                        "4": 1 / 6,
                        "5": 1 / 6,
                        "6": 1 / 6
                    },
                "L":  # 'L' indicates a loaded die
                    {
                        "1": 1 / 10,
                        "2": 1 / 10,
                        "3": 1 / 10,
                        "4": 1 / 10,
                        "5": 1 / 10,
                        "6": 1 / 2
                    }
            }

        self.vtable = [{}]  # we can use this trick to declare a two dimensional list of dicts for the Viterbi table

    def viterbi(self):
        dummy = -999999
        possible_paths = {}

        for state in self.log_states:  # Set states in this initial column to zero probability
            self.vtable[0][state] = dummy

        self.vtable[0]['S'] = 0  # Initialize the start state, as everything must start here, log(1)=0.
        possible_paths['S'] = ['S']

        for position in range(1, self.observed_sequence_length):
            new_paths = {}  # At every new position, replace our previous best path for each state with updated ones
            self.vtable.append({})  # At every new position in the observed sequence, we have a new states dict
            # print(self.vtable)

            for new_state_l in self.log_states:  # iterate over the possible next states

                possible_path_probabilities = []  # We'll make a list of the possibilities, then choose the best one

                for old_state_k in self.log_states:

                    try:
                        # emis = self.emissions[new_state_l][self.observed_sequence[position]] # the 0.3 you see in notes under A
                        # sec = self.states[old_state_k][new_state_l]
                        temp_prob = self.vtable[position - 1][old_state_k] + self.log_states[old_state_k][new_state_l]

                    except KeyError:
                        temp_prob = dummy

                    # self.vtable[position][new_state_l] =
                    # OK THERE IS A BUNCH OF STUFF MISSING HERE, CAN YOU FIGURE IT OUT?
                    # temp_prob = self.vtable[position-1][old_state_k]*sec * emis
                    # print(position, new_state_l, old_state_k)

                    possible_path_probabilities.append((temp_prob, old_state_k))

                    # wrap the state k that gave us this result and the result into a tuple and collect 'em

                # print(new_state_l, ": possible_path_probabilities", sorted(possible_path_probabilities))
                local_max_path_probability, old_state_k_that_gave_local_max_path_probability = max(
                    possible_path_probabilities)

                # print("Local max path probability:", local_max_path_probability, "Old State K:",
                #      old_state_k_that_gave_local_max_path_probability, "New State L", new_state_l)
                # tuples are ordered orthographically, so the second element conveniently "comes along for the ride"
                # now you see why we wrapped up the probability and the state that gave rise to it together in a tuple

                try:
                    emission_prob = self.log_emissions[new_state_l][self.observed_sequence[position]]

                except KeyError:
                    emission_prob = dummy

                self.vtable[position][new_state_l] = local_max_path_probability + emission_prob

                # need to multiply by the emission probability. could have done this prior to maximization,
                # but as it is constant with respect to the iteration this is more efficient.

                new_paths[new_state_l] = possible_paths[old_state_k_that_gave_local_max_path_probability] + \
                                         [new_state_l]
                # print('here', possible_paths[old_state_k_that_gave_local_max_path_probability], new_state_l)

            possible_paths = new_paths
            # print('here', possible_paths)

            # We should now have created a new dict that has as each state as keys, with values
            # that are updated paths for each of those states. Get rid of the out-of-date ones

        return max((self.vtable[self.observed_sequence_length - 1][state],
                    possible_paths[state]) for state in self.states)  # state in self.states

        # List comprehensions are a powerful feature of Python!!
        # This one line accomplishes all of the termination conditions of Viterbi by iterating over the states,
        # and collecting into a list all of the path probabilities found in the last column of the viterbi table
        # actually, it first bundles the probabilities together in a tuple with the state which gave rise to them
        # The max function is used to identify the tuple from the list that has the highest probability.
        # It is this tuple that is returned. The Max function will by default ignore the second element of the tuple
        # unless there are identical first elements (i.e. a tie has occurred). In the case of a tie the tuple containing
        # the state with the lowest alphanumeric sort order should be the one returned.

    def evaluate(self, state_path):
        pro = 0
        for i in range(0, len(state_path) - 1):
            transit_prob = self.log_states[state_path[i]][state_path[i + 1]]
            emission_prob = self.log_emissions[state_path[i]][self.observed_sequence[i]]
            pro = pro + transit_prob
            pro = pro + emission_prob
        pro = pro + (self.log_emissions[state_path[-1]][self.observed_sequence[-1]])
        # you will need to provide some functionality here.
        # The PowerPoint describes the formula and procedure. Can you reduce it to code?
        # Pretty much all of the variables you need are already available to you (at least once you have run Viterbi)
        # ...you have the sequence in self.observed_sequence
        # ...you have the emissions in self.emissions
        # ...you have the transitions in self.states
        return pro


if __name__ == "__main__":
    tuple_acid = FastA_V2.FastA("160_membrane_prots.txt")
    # tuple_acid = FastA_V2.FastA("645_non_membrane_prots.txt.fasta")
    my_states = {
        "S": {
            "+": 0.5,
            "-": 0.5,
        },
        "+": {
            "+": 0.85,
            "-": 0.15,
        },
        "-": {
            "-": 0.95,
            "+": 0.05,
        }
    }

    my_emissions = {

        "S":
            {
                "_": 1
            },
        "-":
            {
                "A": 0.1,
                "C": 0.40,
                "G": 0.40,
                "T": 0.1,

            },
        "+":
            {
                "A": 0.35,
                "C": 0.20,
                "G": 0.10,
                "T": 0.35,
            }
    }

    for annotation, combine in tuple_acid:
        combine_no_space = re.sub(" ", "", combine)
        # print(combine_no_space)
        acid_state = combine_no_space.split("#")
        acid = acid_state[0]
        # try:
        #     acid_state[1]
        # except IndexError:
        #     print('Actual path not given.')

        # try:
        #     acid_state = combine_no_space.split("#")
        #
        #     acid = acid_state[0]
        # path = 'Y'

        # except ValueError:
        #     print('Actual path not given.')
        #     acid = combine_no_space

        # path = 'N'
        acid = "_" + acid
        if "ECOLI" in annotation:  # comment out this line and move the codes below it to correct indent for question 2i
            my_HMM = HMM(acid, my_states, my_emissions)
            probability_of_viterbi_decoded_state_path, viterbi_decoded_state_path = my_HMM.viterbi()
            print(">" + annotation)
            print('Sequence: ', acid)
            fake_path = 'S' + 'I' * (len(acid) - 1)
            # state_path = ''.join(viterbi_decoded_state_path)
            # print("         ", state_path)
            print('Viterbi:  ', ''.join(viterbi_decoded_state_path))
            print("Fake path:", fake_path)
            try:
                print('Actual:   ', 'S' + acid_state[1])
            except IndexError:
                print('Actual path not given.')

            print("Probability Viterbi", probability_of_viterbi_decoded_state_path)
            # probability_eval = my_HMM.evaluate(viterbi_decoded_state_path)
            # print('Evaluate probability:', probability_eval)
            probability_fake_path = my_HMM.evaluate(fake_path)
            print('Probability non-membrane:', probability_fake_path)
            print('Odds ratio:', math.exp(probability_of_viterbi_decoded_state_path - probability_fake_path))
            print("")

# compare evaluate and end max viterbi value from viterbi program
# estimate how good, forward viterbi sum instead of max
# print(math.exp(-1001.2784709812921--1026.6084269854787))
