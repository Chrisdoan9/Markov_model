from pickle import dump
from random import random
from textwrap import fill
from FastA_V2 import FastA
import re

__author__ = 'Wombat'


class Hamlet(object):
    """Hamlet.py  Makes a Markov model from Hamlet's speeches, and outputs new speeches based on that model.
    Version 3.0  (17 April 2020), by Paul Fawcett, but broadly inspired by an earlier PERL program by Jeff Elhai
    Version 3.0 is a conversion from Python 2.x to Python 3.x

    INPUT FILES: Text file containing arbitrary training text, in format:
              {text}
              {text}
        Each brace of text represents one utterance
        In creating new speeches, the program will begin with { and end with } (but won't print either symbol)

    OUTPUT FILES: Hamlet_dict
            a Python pickle file containing a binary serialization of the key transition dict.
    """

    def __init__(self, filepath=None, order=1):

        self.sum_symbol = '|'  # A special symbol used in the keys of the main model dict to indicate a sum for a stem
        start_symbol = 'S'  # All path actually begin with this special symbol
        # end_symbol = '}'        # All speeches actually end with this special symbol

        line = ''  # This will sequentially contain each of the utterances from the input file

        self.order = order  # This corresponds to the order of the Markov chain. See the notes

        self.model = {}

        # The central dict-of-dicts of the program has keys corresponding to words of length order.
        # These keys correspond to stems.  The values are generally also dicts.  Each of these dicts in turn
        # has keys which are also words of length order, but now corresponding to overlapping extensions of the
        # stem. The values in these dicts initially corresponds to a count of the number of times that a particular
        # stem has been observed, but are later converted to frequencies

        self.cumulative_model = {}  # a version of the above dict where the frequencies have been made cumulative

        self.start_string = start_symbol * order  # String that added to the beginning and end of each sequence.
        # self.end_string = end_symbol * order        # These embody the idea that start and end are special states

        if not filepath:  # If no speech file has been specified

            line = 'How much wood could a woodchuck chuck if a woodchuck could chuck wood?'
            line = self.start_string + line  # + self.end_string
            self._determine_raw_counts(line)

        else:  # A speech file has been specified, read it in utterance by utterance

            print("Analysing", filepath)
            sequence_object = FastA(filepath)

            for annotation, line in sequence_object:
                combine_no_space = re.sub(" ", "", line)
                acid_state = combine_no_space.split("#")
                state = acid_state[1]
                if 'ECOLI' not in annotation:
                    line = self.start_string + state  # + self.end_string
                    self._determine_raw_counts(line)

        # Now Convert the raw counts of each transition into frequencies.  As usual, we will be using these frequencies
        # as estimates of the underlying transition probabilities.

        for stem, total_for_stem in self.model.items():  # grab each key / value pair
            # print(stem, total_for_stem)

            if stem[-1] == self.sum_symbol:  # only interested in ones containing a stem count sum!
                predictions = self.model[stem[0:-1]]  # retrieve the corresponding sub-dict with that stem as key
                for prediction, count in predictions.items():  # iterate over the key / value pairs of the sub dict
                    predictions[prediction] = count / total_for_stem  # ..converting the counts in frequencies as we go

        # It will be much easier to work with cumulative distributions when we get into the business of
        # generating speeches, so make a companion dict-of-dict with this structure.
        self.cumulative_model = self._makecumulative(self.model)

    def _determine_raw_counts(self, line):

        """ Takes a line corresponding to an utterance in the input file and uses it to update our dict-of-dicts
         containing the raw counts of how many times one word transitions to another.
        """

        num_words = len(line) - self.order

        for i in range(num_words):

            stem = line[i: i + self.order]
            prediction = line[i + 1: i + self.order + 1]
            # print(i, stem, prediction)

            if stem not in self.model:  # have we seen this stem before? If not:

                self.model[stem] = {prediction: 0}  # create a dict for predictions when a new stem encountered
                self.model[stem + self.sum_symbol] = 0  # create an entry containing the sum for a given stem

            if prediction not in self.model[stem]:  # perhaps we have encountered a stem before, but we have
                # come across a new prediction for that stem
                self.model[stem][prediction] = 0  # counts are set initially to zero, as we always add next

            self.model[stem][prediction] += 1  # Add one regardless of whether this was a new stem or
            self.model[stem + self.sum_symbol] += 1  # prediction, or something we had seen before

    def make_speech(self, max_length=2500, characters_per_line=80):

        if not self.cumulative_model:
            print("You haven't estimated any probabilities yet, so I can't make a speech")
            return

        speech = []
        cur_symbol = self.start_string  # Usually we wish to initialize with our starting symbol
        # but optionally we could just pick a random stem
        for i in range(max_length):
            possibilities = self.cumulative_model[cur_symbol]
            rand_num = random()  # grab a number between 0 and 1

            for possibility, frequency in sorted(possibilities.items()):
                if rand_num < frequency:
                    speech.append(possibility[-1])
                    cur_symbol = possibility
                    break
            # if cur_symbol == self.end_string:
            #     break

        speech = speech[0:-self.order]  # lop off the end-of-sequence string
        speech = ''.join(speech)  # convert from a list to a string

        print(fill(speech, characters_per_line))  # the fill method of textwrap adds EOL characters to a
        # string so that it has no more than some set number of chars per line

    def _makecumulative(self, distribution):  # Utility to convert distribution into cumulative probabilities

        cumulative_distribution = {}
        for i, i_dict in sorted(distribution.items()):
            if i[-1] == self.sum_symbol:
                continue  # ignore (skip over) the odd-ball entries in the dict that correspond to stem sums
            cumulative_prob = 0
            cumulative_distribution[i] = {}
            for j, val in sorted(i_dict.items()):
                cumulative_prob += val
                cumulative_distribution[i][j] = cumulative_prob

        return cumulative_distribution

    def save_model(self, file_path='acid_dict'):  # object serialization with Python's 'Pickle' module
        self.model = {key: val for key, val in self.model.items() if
                      '|' not in key}  # dictionary comprehension with condition
        # for i in list(self.model):
        # del self.model['S|']
        # del self.model['I|']
        # del self.model['M|']
        # del self.model['O|']
        print(self.model)

        with open(file_path, 'wb') as f:
            dump(self.model, f)  # Pickle is very easy to use. Just specify an object to serialize, and an open
            # file handle, and you are ready to rock.
        print("We have pickled the dict to file", file_path)
        print(self.model)


def main():
    my_speech_maker = Hamlet("160_membrane_prots.txt", 1)  # I should change the name of class but I think you still understand.
    make_another = True
    while make_another:
        my_speech_maker.make_speech()
        user_input = input('\nShall I make another path: ')
        print("\n")
        user_input = user_input.upper()

        try:
            if user_input[0] != 'Y':
                make_another = False
        except IndexError:  # This catches the case where the user just hits return with no other input
            make_another = False

    my_speech_maker.save_model()


if __name__ == '__main__':
    main()
