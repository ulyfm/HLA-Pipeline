import numpy as np
import pandas as pd

aa_dict = {"A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5,
           "H": 6, "I": 7, "K": 8, "L": 9, "M": 10, "N": 11,
           "P": 12, "Q": 13, "R": 14, "S": 15, "T": 16, "V": 17,
           "W": 18, "Y": 19}


def convert_aa(letter: str) -> int:
    """
    Converts a letter to an index based on amino acids as an alphabet
    TODO: error handling
    """
    return aa_dict[letter]


class DefragSTree:
    """
    Stores peptide data combined as a suffix tree in order to efficiently determine whether a given peptide is a
    fragment. Addition of peptide data is lazy, in order to allow the user to add it sequentially (i.e. first
    very large sequences and then shorter ones) so that the same tree can be reused.

    Retention time (RT) is stored at every level of the tree at least as long as min_length for every suffix or prefix
    of a suffix.
    """

    def __init__(self, items: pd.DataFrame, min_length=6, max_length=22):
        """
        Constructor.
        :param items: A dataframe with all relevant peptides, including HLAP_sequence, HLAP_length, HLAP_RT, and
                      HLAP_fragment (initialized to False is ok).
        :param min_length: The minimum length of peptide to consider for fragment removal.
        :param max_length: The maximum length of peptide to consider for fragment removal. (Even longer peptides may
                           still be used to determine whether shorter ones are fragments.)
        """
        self.tree = ([], [None] * 20)  # Tuple of (suffix tree list, RT list)
        self.data = items
        self.min_length = min_length
        self.max_length = max_length
        self._curr_min_length = 1000  # TODO: magic number
        self._length = 0

    def add_above_length(self, length: int):
        """
        Adds all peptides above the given length (and their suffixes) to the suffix tree.
        :param length: length to add peptides ABOVE (non inclusive)
        """
        filtered = self.data[(self.data['HLAP_length'] > length) & (self.data['HLAP_length'] <= self._curr_min_length)]
        items = filtered['HLAP_sequence']
        rt_values = filtered['HLAP_RT']
        for item, rt in zip(items.tolist(), rt_values.tolist()):
            self.add_with_suffixes(item, rt)
            self._length += 1
        self._curr_min_length = length

    def add_with_suffixes(self, item: str, rt: float):
        """
        Adds a peptide and its RT values to the suffix tree, including all its suffixes longer than min_length.
        """
        for i in range(len(item) - self.min_length + 1):
            self.add(item[i:], rt)

    def add(self, item: str, rt: float):
        """
        Adds a string and its RT value to the suffix tree, with the RT value being stored at every level longer than
        min_length.
        """
        depth = 0
        level = self.tree
        for char in item:
            index = convert_aa(char)
            depth += 1
            if level[1][index] is None:
                level[1][index] = ([rt] if depth >= self.min_length else [], [None] * 20)
            else:
                if depth >= self.min_length:
                    level[1][index][0].append(rt)
            level = level[1][index]

    def get(self, key: str):
        """
        Efficiently determines if key is a substring of any of the peptides already added to the tree,
        and returns all the corresponding RT values if so (otherwise None)
        :param key: the string to search for
        """
        level = self.tree
        for char in key:
            index = convert_aa(char)
            if level[1][index] is None:
                return None
            level = level[1][index]
        # print("For", key, "returning", level[0])
        return level[0]


"""
This code (defrag and find_fragments) is partially based on code privately provided by Andrew Fiore-Gartland.
The format of the code is similar but the algorithm has been changed significantly (to use suffix
trees for efficiency).
"""


def defrag(d: pd.DataFrame, rt_thresh=0.5, max_len=22, min_len=6):
    """
    Determines whether peptides are fragments. Peptides from min_len to max_len are processed. A peptide is considered
    a fragment if it is a substring of a longer fragment with an RT difference less than rt_thresh.
    :param d: dataframe to modify
    :param rt_thresh: the maximum absolute difference in retention time (RT) for a peptide to be a fragment
    :param max_len: the max length of peptide to defrag (even longer lengths are compared to but not defragged)
    :param min_len: the min length of peptide to defrag
    """
    d = d.copy()
    d.loc[:, 'HLAP_fragment'] = False
    stree = DefragSTree(d, min_length=min_len, max_length=max_len)
    for i in range(max_len, min_len - 1, -1):
        find_fragments(d, stree, short_len=i, rt_thresh=rt_thresh)
    return d


def find_fragments(d, stree: DefragSTree, short_len=8, rt_thresh=0.5):
    """
    Finds and updates the HLAP_fragment column for peptides of a particular length (short_len)
    :param d: table of peptides and relevant columns
    :param stree: suffix tree to update and search for longer peptides. n.b. it is permanently modified by this method.
    :param short_len: length of peptide for which to look for fragments
    :param rt_thresh: maximum absolute RT diff to allow for a peptide that is a substring of a longer one
    """
    stree.add_above_length(short_len)

    def _check(r):
        seq_rts = stree.get(r['HLAP_sequence'])
        if seq_rts is not None and any(np.abs(rt - r['HLAP_RT']) < rt_thresh for rt in seq_rts):
            return True
        return False

    filtered_d = d[d['HLAP_length'] == short_len]
    d.loc[filtered_d.index, 'HLAP_fragment'] = d.loc[filtered_d.index, 'HLAP_fragment'] | filtered_d.apply(_check, axis=1)

