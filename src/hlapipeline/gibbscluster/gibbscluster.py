# Alphabet and uniprot frequency copied from compute_logo.py in palmotif by Andrew Fiore-Gartland.
import math
import random
from copy import deepcopy
from collections import OrderedDict
import blosum

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

def _normalized_reference_frequencies(reference_frequencies):
    normalized_reference_frequencies = {}
    total = 0
    for aa in reference_frequencies.keys():
        total += reference_frequencies[aa]
    for aa in reference_frequencies.keys():
        normalized_reference_frequencies[aa] = reference_frequencies[aa] / total
    return normalized_reference_frequencies

normalized_uniprot = _normalized_reference_frequencies(uniprot_frequency)

# Input: peptides list with length >= n
# Input: target length n
# Input: reference frequencies (optional)
# Output: array of original strings as tuples along with the offset corresponding to the alignment and a copy of the
#         offset for uses where the previous state of the list is stored.
def gibbs_cluster(peptides: list[str], n: int, reference_frequencies=None) -> list[list[str, int, int]]:
    print(peptides)
    if reference_frequencies is None:
        reference_frequencies = uniprot_frequency
    normalized_reference_frequencies = _normalized_reference_frequencies(reference_frequencies)
    # assemble list of lists
    alignment_list = []
    for peptide in peptides:
        if len(peptide) >= n:
            random_start_index = random.randint(0, len(peptide) - n)
            alignment_list.append(
                [peptide, random_start_index, random_start_index])  # want to do this multiple times?
        else:
            print("skipping short sequence", peptide, "in logo generation")
    # Now that the list
    current_E = _energy_of_alignment(alignment_list, 9, normalized_reference_frequencies)
    moves = 5000
    T = 0.15
    T_end = 0.01
    T_steps = 10
    # Generate a list of steps so that start and end can be defined, due to division inaccuracies
    steps_list = [T]
    step_size = (T - T_end) / (T_steps - 1)
    for i in range(1, T_steps - 1):
        steps_list.append(T - step_size * i)
    steps_list.append(T_end)
    print(steps_list)

    for step in steps_list:
        print("starting at T=", step)
        for j in range(moves):
            if random.random() >= 0.001:  # move type 1 (single sequence move) - 1 in 1000
                move_index = random.randrange(0, len(alignment_list))
                alignment_list[move_index][2] = alignment_list[move_index][1] # Save old value
                alignment_list[move_index][1] = random.randint(0, len(alignment_list[move_index][0]) - n)
                new_E = _energy_of_alignment(alignment_list,  9, normalized_reference_frequencies)
                try:
                    P = min(1.0, math.exp((new_E - current_E) / step))
                except OverflowError:
                    print("newE", new_E, "currentE", current_E, "step", step)
                    P = 0
                if random.random() <= P:
                    current_E = new_E
                else:
                    alignment_list[move_index][1] = alignment_list[move_index][2]  # Restore old value
            else:  # move type 2 (frame shift move)
                move_amount = random.randint(0, 18) - 9
                for i in range(len(alignment_list)):
                    # print(new_alignment_list[i])
                    if move_amount < 0:
                        most_negative = -alignment_list[i][1]
                        negative = max(most_negative, move_amount)
                        alignment_list[i][2] = alignment_list[i][1]
                        alignment_list[i][1] += negative
                    elif move_amount > 0:
                        most_positive = len(alignment_list[i][0]) - n - alignment_list[i][1]
                        positive = min(most_positive, move_amount)
                        alignment_list[i][2] = alignment_list[i][1]
                        alignment_list[i][1] += positive

                new_E = _energy_of_alignment(alignment_list, 9, normalized_reference_frequencies)
                P = min(1.0, math.exp((new_E - current_E) / step))
                if random.random() <= P:
                    current_E = new_E
                else:
                    for i in range(len(alignment_list)):
                        alignment_list[i][1] = alignment_list[i][2]

    print(alignment_list)
    return alignment_list

blosum62 = blosum.BLOSUM(62)
blosum62_probability = {}
_aa_table = dict.fromkeys(uniprot_frequency, 0)
sum_ab = 0
for _aa in _aa_table.keys():
    blosum62_probability[_aa] = dict()
    for _bb in _aa_table.keys():
        # print(aa, bb)
        _Sab = blosum62[_aa][_bb]
        _Qb = normalized_uniprot[_bb]
        _Qa = normalized_uniprot[_aa]
        _Qab = _Qb * _Qa * math.exp(_Sab * 0.31723)  #
        blosum62_probability[_aa][_bb] = _Qab
        print(_Qab)
        sum_ab += _Qab
        # print("qab=", Qab)
print(sum_ab)


def _energy_of_alignment(alignment_list: list[list[str, int, int]], n: int,
                         reference_frequencies: dict) -> float:
    # Weighting: first calculate the number of different amino acids at all positions:
    pos_count = []
    for i in range(n):
        pos_count.append({})
    for i in range(n):
        for [sequence, index, old_index] in alignment_list:
            res = sequence[i + index]
            if res not in pos_count[i]:
                pos_count[i][res] = 1
            else:
                pos_count[i][res] += 1
    alpha = 20
    beta = 50  # magic number for Henikoff method
    E = 0
    for i in range(n):
        aa_table = dict.fromkeys(reference_frequencies, 0)  # use custom reference?
        total_aas = 0
        for j in range(len(alignment_list)):
            [sequence, index, old_index] = alignment_list[j]
            aa_table[sequence[index + i]] += 1
            total_aas += 1
        for aa in aa_table.keys():
            Cpa = aa_table[aa]
            Fpa = aa_table[aa] / total_aas
            Gpa = 0
            for bb in aa_table.keys():
                Fpb = aa_table[bb] / total_aas
                Qab = blosum62_probability[aa][bb]
                Qb = reference_frequencies[bb]
                Gpa += Fpb / Qb * Qab
            Ppa = (alpha * Fpa + beta * Gpa) / (alpha + beta)
            Qa = reference_frequencies[aa]
            E += Cpa * math.log(Ppa / Qa)  # correct base?
    return E


# assuming sequence list does not contain duplicates.
def _generate_weights(sequence_list: list[str], threshold: float) -> dict[str, float]:
    pass

