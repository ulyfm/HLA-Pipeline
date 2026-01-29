import os
import re

import matplotlib.pyplot as plt
import seaborn as sns
from palmotif import compute_motif, svg_logo, uniprot_frequency

from hlapipeline import LocalDatabase
from hlapipeline.defrag import defrag
from hlapipeline.files.OutputAdapter import OutputAdapter
from hlapipeline.files.PeptideFileAdapter import PeptideFileAdapter


class HLAPipeline:
    """
    Overall HLA pipeline. Manages data processing steps in order.
    """

    def __init__(self, output: OutputAdapter, peptide_file: PeptideFileAdapter, local_db: LocalDatabase = None) -> None:
        """
        Constructor. Loads data from the peptide file and records key metadata.
        :param output: the output adapter (where files will be output to)
        :param peptide_file: the input file adapter (where data will be pulled from)
        :param local_db: the local database with additional gene information to optionally join to the main data
        """
        self._db_access = local_db
        self.final_row = {}
        self._output = output
        self._peptide_file = peptide_file
        self._data_table = peptide_file.get_dataframe()
        self.final_row['file_name'] = peptide_file.get_file_name()
        self._name_before_group = peptide_file.get_base_name()
        self.final_row['date_created'] = peptide_file.get_date_created()
        self.final_row['HLA_allele'] = peptide_file.get_allele()

    def remove_sp(self) -> None:
        """
        Remove and save peptides with 'sp' in the 'Master Protein Accessions' column. 'sp' peptides are considered
        junk peptides.
        Removed peptides are saved in a folder called sp_peptides.

        :return: the starting table, minus peptides with sp.
        """
        self.final_row['sp_count'] = sum(self._data_table['HLAP_master_accessions'].str.contains("sp", na=False))
        print("Removing", self.final_row['sp_count'], "sp peptides")
        self._output.save_csv_infer_name(self._data_table[self._data_table['HLAP_master_accessions']
                                         .str.contains("sp", na=False)], self._name_before_group, "sp_peptides")
        self._data_table = self._data_table[~self._data_table['HLAP_master_accessions'].str.contains("sp", na=False)]
        self._output.save_csv_infer_name(self._data_table, self._name_before_group, "spRM_peptides")

    def remove_frag(self) -> None:
        """
        Uses a special defrag script to analyze whether peptides are fragments. The results are stored in a new column
        and then fragments are dropped from the table and stored in a backup file.

        A peptide is a fragment if there is a longer peptide that it is a substring of AND the difference in RT is
        less than a certain threshold (default 0.5)

        :return: the starting dataframe, minus rows that are frags and with an additional boolean column 'fragment'
        """
        if 'HLAP_RT' not in self._data_table.columns:
            print("Skipping frag removal due to no RT column.")
            return
        nosp_table = defrag.defrag(self._data_table)
        self.final_row['fragment_count'] = sum(nosp_table['HLAP_fragment'])  # Final
        self._output.save_csv_infer_name(nosp_table[nosp_table['HLAP_fragment']], self._name_before_group,
                                         "frag_peptides")
        self._data_table = nosp_table[~nosp_table['HLAP_fragment']]
        self._output.save_csv_infer_name(self._data_table, self._name_before_group, "fragRM_peptides")
        print("Removing", self.final_row['fragment_count'], "fragment peptides")

    def remove_duplicates(self):
        """
        Updates the dataframe to not have any duplicates (keeping only one). Duplicates removed are stored in a backup
        file.
        """
        duplicates = self._data_table.duplicated(subset='HLAP_sequence', keep="last")
        dup_table = self._data_table[duplicates]
        self._output.save_csv_infer_name(dup_table, self._name_before_group, "dup_peptides")
        nodup_table = self._data_table[~duplicates]
        self.final_row['duplicate_count'] = self._data_table.shape[0] - nodup_table.shape[0]
        print("Removing", self.final_row['duplicate_count'], "peptides with duplicate sequences")
        self.final_row['total_peptides'] = nodup_table.shape[0]
        print("There are now", self.final_row['total_peptides'], "peptides remaining")
        self._data_table = nodup_table

    def _generate_mers_piechart(self):
        """
        Generates two pie charts for a given table, based on peptide length. One pie chart includes only 8-14mers.
        The other pie chart has 7-14mers as well as an 'other' category.
        """
        mers_labels = []
        mers_cts = []
        total_mers = 0
        for i in range(8, 15):
            mers_labels.append(str(i) + "mers")
            sum_mer = sum(self._data_table['HLAP_length'] == i)
            total_mers += sum_mer
            mers_cts.append(sum_mer)
            self.final_row[str(i) + "mers"] = sum(self._data_table['HLAP_length'] == i)
        if set(mers_cts) == {0}:
            return
        fig, ax = plt.subplots()
        ax.pie(mers_cts, labels=mers_labels, colors=sns.color_palette("Set3"), autopct='%1.1f%%')
        self._output.save_fig_infer_name(fig, self._name_before_group, "FragRM_spRM_PeptideGroups_PieChart_Abridged")
        plt.close(fig)
        mers_labels.append("other")
        mers_cts.append((self.final_row['total_peptides'] - total_mers))
        mers_labels.append("7mers")
        mers_cts.append(sum(self._data_table['HLAP_length'] == 7))
        fig, ax = plt.subplots()
        ax.pie(mers_cts, labels=mers_labels, colors=sns.color_palette("Set3"), autopct='%1.1f%%')
        self._output.save_fig_infer_name(fig, self._name_before_group, "FragRM_spRM_PeptideGroups_PieChart_FULL")
        plt.close(fig)

    def identify_cotransduced(self, peptide: str = ''):
        """
        Queries the user for a search term for con-transduced peptides and then saves the search results to the final
        table.
        :param peptide: an optional pre-entered peptide to use instead of user input.
        """
        if peptide == '':
            print("Co-transduced peptides - comma separated? (enter to skip)")
            peptide = input("> ")
        if peptide != '' and peptide != 'none':
            if "," in peptide:
                peptide_spl = peptide.split(",")
                num = 0
                for i in range(0, len(peptide_spl)):
                    num = self._process_cotrans_peptide(peptide_spl[i].strip(), num)
            else:
                self._process_cotrans_peptide(peptide, 0)

    def _process_cotrans_peptide(self, peptide: str, num: int):
        """
        Queries the table for peptides with matching accessions or descriptions to the "co-transduced peptide"
        specified. Matching peptides are added as a list to the final table.

        :param peptide: the string peptide to search for (using regular expressions, with caveats -- see README)
        :param num: the number to list the peptide as and add to the final table, if multiple searches are being done.
        """

        def _clean_peptide(dirty_peptide: str, match_pept: str):
            dirty_peptide = dirty_peptide.lower()
            if "-" not in match_pept:
                dirty_peptide = dirty_peptide.replace("-", "")
            if "_" not in match_pept:
                dirty_peptide = dirty_peptide.replace("_", "")
            if " " not in match_pept:
                dirty_peptide = dirty_peptide.replace(" ", "")
            return dirty_peptide.strip()

        def _match_peptide(data_pept: str, match_pept: str) -> bool:
            if data_pept is None or match_pept is None:
                return False
            match_pept = match_pept.lower()
            data_pep_spl = []
            try:
                data_pep_spl = data_pept.split(";")
            except AttributeError:
                print("Data is invalid, skipping match row", data_pept)
            for i in range(0, len(data_pep_spl)):
                if bool(re.fullmatch(match_pept, _clean_peptide(data_pep_spl[i], match_pept))):
                    return True
            return False

        peptide = _clean_peptide(peptide, peptide).lower()
        bool_table = self._data_table['HLAP_master_accessions'].apply(_match_peptide, match_pept=peptide)
        bool_table |= self._data_table['HLAP_master_descriptions'].apply(_match_peptide, match_pept=peptide)
        bool_table |= self._data_table['HLAP_protein_accessions'].apply(_match_peptide, match_pept=peptide)
        peptide1_matches = self._data_table[bool_table]
        list_match = peptide1_matches[
            ['HLAP_sequence', 'HLAP_master_accessions', 'HLAP_master_descriptions',
             'HLAP_protein_accessions']].values.tolist()
        match_dict = {}
        for row in list_match:
            for i in range(1, 4):
                data_pep_spl = row[i].split(";")
                for pep in data_pep_spl:
                    pep = pep.strip()
                    if bool(re.fullmatch(peptide, _clean_peptide(pep, peptide))):
                        if pep not in match_dict:
                            match_dict[pep] = set()
                        match_dict[pep].add(row[0])
                        # break - trying to match all...
        match_count = num
        peptide_count = 0
        for key in match_dict.keys():
            self.final_row['co-transduced_protein_' + str(match_count)] = key
            self.final_row['co-transduced_peptides_' + str(match_count)] = ", ".join(sorted(list(match_dict[key])))
            # This sorted() in the line above makes the output deterministic, which helps with testing and comparisons
            peptide_count += len(match_dict[key])
            match_count += 1
        print("-", peptide, "matches:", len(peptide1_matches), "peptides, or", peptide_count,
              "including multiple matches; in", match_count - num, "groups")
        return match_count

    def _generate_logos(self):
        """
        Generates logo images (showing common amino acids at different positions in the peptides, scaled by frequency).
        Images are exported to <base directory>/logo_results/<HLA file name>/<n>mer_logo.png
        """
        for i in range(8, 15):
            mers = self._data_table[self._data_table['HLAP_length'] == i]
            seqs = mers['HLAP_sequence'].tolist()
            if not seqs:
                continue
            # print(seqs, ",", uniprot_frequency)
            motif = compute_motif(seqs, reference_freqs=uniprot_frequency) \
                .rename(columns=dict(zip(list(range(0, i)), list(range(1, i + 1)))))
            logo_dir = self._output.init_logo_directory(self._name_before_group)
            svg_logo(motif, os.path.join(logo_dir, str(i) + "mer_logo.svg"), color_scheme='chemistry')

    def _db_data(self):
        """
        Appends ProteinAtlas / Gene Ontology data to the pipeline table by left joining on accession codes.
        For accession codes with multiple elements separated by semicolons, additional columns will be generated
        up to the specified depth. For example, if there are three accession codes (e.g. "ABCD; EFGH; IJKL") and a
        depth of 2, the data table will have the PA/GO data for ABCD and EFGH, appended as columns with suffixes _0
        and _1, in order.
        """
        if not self._db_access:
            print("Skipping 3rd party database addition features.")
            return
        self._data_table['HLAP_accession'] = self._data_table['HLAP_master_accessions']
        self._data_table = self._data_table.fillna('')

        def _split_accession(acc_list: str, index: int):
            """
            Utility method to split accessions. The accession code at the given index is returned.
            (for future: spaces can be problematic -- could also do .strip instead of splitting by "; ")
            :param acc_list: string list of accession codes to split.
            :param index: the index of the accession code to pull from the string.
            :return: the accession code at the given index.
            """
            spl = acc_list.split("; ")
            if index < len(spl):
                return spl[index]
            else:
                return ''

        for i in range(0, self._data_table['HLAP_accession'].str.split("; ")):
            self._data_table['HLAP_accession' + str(i)] = self._data_table['HLAP_accession'].apply(
                lambda x: _split_accession(x, i))
            self._data_table = self._data_table.merge(self._db_access.get_table(), left_on='HLAP_accession' + str(i),
                                                      right_on='HLAP_accession', how='left',
                                                      suffixes=(None, "_" + str(i)))
            self._data_table.drop('HLAP_accession_' + str(i), axis=1, inplace=True)

    def get_result(self, skip_cleanup: bool, skip_cotransduced: bool, assumecotransduced: bool, skip_dataviz: bool,
                   cotransduced_peptide='', database_features=False) -> dict:
        """
        Retrieve the final metadata row. This includes things such as name, date, the numbers of peptides excluded,
        and the numbers of 7-14mers included.

        :return: the final metadata row (dict) as constructed by this pipeline.
        """
        print("File contains", len(self._data_table), "peptides")
        if not skip_cleanup:
            self.remove_sp()
            self.remove_frag()
            self.remove_duplicates()
        if assumecotransduced:
            self._process_cotrans_peptide(self._peptide_file.get_cotransduced(), 0)
        self._output.save_csv_infer_name(self._data_table, self._name_before_group, "dupRM_fragRM_spRM_peptides")
        if database_features:
            self._db_data()
        # Double check AGAIN that there aren't duplicates
        assert len(self._data_table) == len(self._data_table.drop_duplicates(subset="HLAP_sequence"))
        self._output.save_csv_infer_name(self._data_table, self._name_before_group, "final_peptides")
        self._output.save_csv_infer_name(self._data_table[(self._data_table['HLAP_length'] >= 8)
                                         & (self._data_table['HLAP_length'] <= 14)],
                                         self._name_before_group, "final_peptides_8-14")
        if not skip_cotransduced:
            self.identify_cotransduced(cotransduced_peptide)
        if not skip_dataviz:
            self._generate_mers_piechart()
            self._generate_logos()
        return self.final_row

    def get_data_table(self):
        """
        Returns CLONE of data table in progress.
        :return:
        """
        return self._data_table.copy()
