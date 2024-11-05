# Alphabet and uniprot frequency copied from compute_logo.py in palmotif by Andrew Fiore-Gartland.
import math
import random
from copy import deepcopy
from collections import OrderedDict

aa_alphabet = [aa for aa in 'ARNDCQEGHILKMFPSTWYVBZX-']

uniprot_frequency = {'A': 8.25,
                     'R': 5.53,
                     'N': 4.06,
                     'D': 5.45,
                     'C': 1.37,
                     'Q': 3.93,
                     'E': 6.75,
                     'G': 7.07,
                     'H': 2.27,
                     'I': 5.96,
                     'L': 9.66,
                     'K': 5.84,
                     'M': 2.42,
                     'F': 3.86,
                     'P': 4.70,
                     'S': 6.56,
                     'T': 5.34,
                     'W': 1.08,
                     'Y': 2.92,
                     'V': 6.87}


# Input: peptides list with length >= n
# Input: target length n
# Input: reference frequencies (optional)
# Output: array of original strings as tuples along with the offset corresponding to the alignment
def gibbs_cluster(peptides: list[str], n: int, reference_frequencies=None) -> list[tuple[str, int]]:
    print(peptides)
    if reference_frequencies is None:
        reference_frequencies = uniprot_frequency
    normalized_reference_frequencies = {}
    total = 0
    for aa in reference_frequencies.keys():
        total += reference_frequencies[aa]
    for aa in reference_frequencies.keys():
        normalized_reference_frequencies[aa] = reference_frequencies[aa] / total
    print(normalized_reference_frequencies)
    # assemble list of tuples
    alignment_list = []
    for peptide in peptides:
        if len(peptide) >= n:
            alignment_list.append(
                (peptide, random.randrange(0, len(peptide) - n + 1)))  # want to do this multiple times?
        else:
            print("skipping short sequence", peptide, "in logo generation")
    alignment_index = 0
    current_E = _energy_of_alignment(alignment_list, alignment_index, 9, normalized_reference_frequencies)
    moves = 5000
    T = 0.15
    T_end = 0.01
    T_curr = T
    T_steps = 10

    while T_curr >= T_end:
        print("starting at T=", T_curr)
        for j in range(moves):
            if random.random() >= 0.00:  # move type 1 (single sequence move)
                move_index = random.randrange(0, len(alignment_list))
                prev_val = alignment_list[move_index]
                move_frame_index = random.randrange(0, len(prev_val[0]) - n + 1)
                alignment_list[move_index] = (prev_val[0], move_frame_index)
                new_E = _energy_of_alignment(alignment_list, 0, 9, normalized_reference_frequencies)
                P = min(1.0, math.exp((new_E - current_E) / T_curr))
                if random.random() <= P:
                    current_E = new_E
                else:
                    alignment_list[move_index] = prev_val

            else:  # move type 2 (frame shift move)
                move_amount = random.randint(0, 18) - 9
                # print("Type 2 rand", move_amount)
                new_alignment_list = deepcopy(alignment_list)
                for i in range(len(new_alignment_list)):
                    # print(new_alignment_list[i])
                    if move_amount < 0:
                        most_negative = -new_alignment_list[i][1]
                        negative = max(most_negative, move_amount)
                        # print("negative", negative, "most negative", most_negative)
                        # print("new loc", new_alignment_list[i][1] + negative, "seq length", len(new_alignment_list[i][0]))
                        new_alignment_list[i] = (new_alignment_list[i][0], new_alignment_list[i][1] + negative)
                    elif move_amount > 0:
                        most_positive = len(new_alignment_list[i][0]) - n - new_alignment_list[i][1]
                        positive = min(most_positive, move_amount)
                        # print("positive", positive, "most positive", most_positive)
                        # print("new loc", new_alignment_list[i][1] + positive, "seq length",
                        #       len(new_alignment_list[i][0]))
                        new_alignment_list[i] = (new_alignment_list[i][0], new_alignment_list[i][1] + positive)

                new_E = _energy_of_alignment(new_alignment_list, 0, 9, normalized_reference_frequencies)

                if new_E != current_E:
                    P = min(1.0, math.exp((new_E - current_E) / T_curr))
                    if random.random() <= P:
                        print("Doing move type 2")
                        current_E = new_E
                        alignment_list = new_alignment_list
                    else:
                        # print("Not doing move type 2")
                        pass
                    print("Move type 2 prob", P, "from newE", new_E, "and oldE", current_E)

        T_curr -= (T - T_end) / 10
    print(alignment_list)
    return alignment_list


def _energy_of_alignment(alignment_list: list[tuple[str, int]], alignment_index: int, n: int,
                         reference_frequencies: dict) -> float:
    E = 0
    for i in range(n):
        aa_table = dict.fromkeys(reference_frequencies, 0)  # use custom reference?
        for (sequence, index) in alignment_list:
            # print(sequence, index, alignment_index, i)
            loc = index + alignment_index + i
            # print("Sequence length", len(sequence), "loc=", loc, "seq=", sequence)
            aa_table[sequence[loc]] += 1
        for aa in aa_table.keys():
            Cpa = aa_table[aa]
            Ppa = (aa_table[aa] + 1) / len(
                alignment_list)  # TODO - no sequence weighting or pseudocount correction for now; just adding 1
            Qa = reference_frequencies[aa]
            E += Cpa * math.log(Ppa / Qa)  # correct base?
    return E


# assuming sequence list does not contain duplicates.
def _generate_weights(sequence_list: list[str], threshold: float) -> dict[str, float]:

    sequence_list.sort()
    clusters = OrderedDict()
   


