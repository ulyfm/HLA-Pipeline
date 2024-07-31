import os.path

import pandas as pd

from hlapipeline.files.PeptideFileAdapterException import PeptideFileException


class PeptideFileAdapter:
    """
    Provides an extensible interface that can be used to standardize various different peptide input files.
    For example, if a different version of ProteomeDiscoverer is used, you can write a subclass that fixes it.
    There are several columns in particular that must be included. See the get_dataframe method for more info.
    """

    def __init__(self, absolute_file_path, allele: str = None):
        """
        Constructor. Sets important values for file adapter
        """
        self._data = None
        self._path = absolute_file_path
        self._allele = allele

    def get_dataframe(self) -> pd.DataFrame:
        """
        Standardizes the dataframe by loading it as a TSV and then adjusting/renaming columns that are required for
        further processing.
        The required columns are (in the format 'original in my format -> required standardized name'):
        Sequence -> HLAP_sequence
        RT -> HLAP_RT
        Master Protein Accessions -> HLAP_master_accessions
        Master Protein Descriptions -> HLAP_master_descriptions
        Length -> HLAP_length
        :throws PeptideFileException if the unicode formatting is not UTF-8 or UTF-16.
        """
        if not os.path.exists(self._path):
            raise PeptideFileException("File path " + str(self._path) + " does not exist.")
        try:
            self._data = pd.read_csv(self._path, delimiter="\t")
        except UnicodeDecodeError:
            try:
                self._data = pd.read_csv(self._path, delimiter="\t", encoding='UTF-16')
            except UnicodeDecodeError:
                raise PeptideFileException("Could not process " + str(self._path) + " because the encoding is "
                                                                                    "unreadable. Only UTF-8 and "
                                                                                    "UTF-16 are accepted.")

        def _clean_seq(an_seq):
            if an_seq == '':
                return ''
            part1 = an_seq[an_seq.index(".") + 1:]
            part2 = part1[:part1.index(".")]
            return part2

        # Sets HLAP_sequence from 'Sequence' or alternatively 'Annotated Sequence'. This is a REQUIRED column.
        if 'Sequence' in self._data.columns:
            self._data['HLAP_sequence'] = self._data['Sequence'].astype(str)
        elif 'Annotated Sequence' in self._data.columns:
            self._data['HLAP_sequence'] = self._data['Annotated Sequence'].astype(str).apply(_clean_seq)
        else:
            raise PeptideFileException("Could not process " + str(self._path) + "because a sequence column could not "
                                                                                "be found.")

        # Attempt to find RT column. If absent, defragmentation cannot occur.
        for name in self._data.columns:
            if name.startswith("RT"):
                self._data['HLAP_RT'] = self._data[name]
                break
        if 'HLAP_RT' not in self._data.columns:
            print("Warning: " + str(self._path) + " is missing an RT column.")

        # Attempt to add accession columns
        if 'Master Protein Accessions' in self._data.columns:
            self._data['HLAP_master_accessions'] = self._data['Master Protein Accessions'].astype(str)
        else:
            self._data['HLAP_master_accessions'] = ""
            print("Warning: " + str(self._path) + " is missing a master accessions column.")

        if 'Protein Accessions' in self._data.columns:
            self._data['HLAP_protein_accessions'] = self._data['Protein Accessions'].astype(str)
        else:
            self._data['HLAP_protein_accessions'] = ""
            print("Warning: " + str(self._path) + " is missing a protein accessions column.")

        # Attempt to add description column.
        if 'Master Protein Descriptions' in self._data.columns:
            self._data['HLAP_master_descriptions'] = self._data['Master Protein Descriptions'].astype(str)
        else:
            self._data['HLAP_master_descriptions'] = ""
            print("Warning: " + str(self._path) + " is missing a master descriptions column.")

        # Add a length column, either one that is already there or by processing HLAP_sequence.
        if 'length' in self._data.columns:
            self._data['HLAP_length'] = self._data['length'].astype(int)
        else:
            self._data = self._data.assign(HLAP_length=self._data['HLAP_sequence'].map(len))
        return self._data

    def get_allele(self) -> str:
        """
        Returns the HLA allele associated with the file. In my case, it can usually be inferred from the file name.
        """
        if self._allele is not None:
            return self._allele
        name_spl = self.get_base_name().split("_")
        if name_spl[2] == "HLA":
            return "HLA_" + str(name_spl[3])
        else:
            print("Could not infer HLA allele for " + self._path)
            allele = input("please specify > ")
            return allele.strip()

    def get_file_name(self) -> str:
        """
        Returns the filename overall (including extension)
        """
        return os.path.basename(self._path)  # Final

    def get_base_name(self) -> str:
        """
        Returns the 'base name', i.e. the file name without unnecessary extensions at the end.
        """
        filename = self.get_file_name()
        if filename.endswith("_PeptideGroups.txt"):
            return filename[:-len("_PeptideGroups.txt")]  # only if peptide group file...
        elif filename.endswith(".txt"):
            return filename[:-len(".txt")]
        else:
            return filename

    def get_date_created(self) -> str:
        """
        Returns the data the file was created, in my case inferred from the file name.
        """
        return self.get_base_name().split("_")[1]

    def get_cotransduced(self) -> str:
        """
        Returns the cotransduced peptides that can be inferred from the file name, IF it follows the correct format.
        """
        spl = self.get_base_name().split("_")
        if spl[4] != 'bRP' and 'bRP' in spl:
            return '-'.join(spl[4:spl.index('bRP')])
        else:
            return ''
