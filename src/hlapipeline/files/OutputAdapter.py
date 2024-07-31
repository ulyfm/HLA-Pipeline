import os

import pandas as pd
from matplotlib.figure import Figure


class OutputAdapter:
    """
    Manages output file access for HLA-Pipeline. The relevant files are the overall output directory, the union table,
    and the overview table.
    """
    def __init__(self, output_dir: str, union_table: str, overview_table: str):
        self._output_dir = output_dir
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        # Overview table
        self._final_table_path = overview_table
        if os.path.exists(overview_table):
            print("Previous final table found, going to append new data...")
            self._final_table = pd.read_csv(overview_table, index_col=0)
        else:
            self._final_table = pd.DataFrame()
        # Union file
        self._union_table_path = union_table
        if os.path.exists(union_table):
            self._union_table = pd.read_csv(union_table, index_col=0)  # TODO: harmonize data types
            if 'sequence' in self._union_table: # Old version backwards compatibility
                print("Renaming union table columns to match current version")
                self._union_table.rename(columns={"sequence": "HLAP_sequence", "length": "HLAP_length",
                                                  "Master Protein Descriptions": "HLAP_master_descriptions"},
                                         inplace=True)
        else:
            self._union_table = pd.DataFrame()

    def save_csv_infer_name(self, dataframe: pd.DataFrame, file_base_name: str, directory_name: str):
        """
        Saves a CSV to the given directory with an inferred file name that looks like file_base_name_directory_name.csv
        :param dataframe: data to save
        :param file_base_name: base file name, e.g. Test_010101_HLA_A_Cov_bRP
        :param directory_name: directory name to save the file to AND append to the file name
        """
        dir_path = os.path.join(self._output_dir, directory_name)
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)
        dataframe.to_csv(os.path.join(dir_path, file_base_name + "_" + directory_name + ".csv"))

    def save_fig_infer_name(self, fig: Figure, file_base_name: str, figname: str):
        """
        Saves a CSV to the image_output directory with an inferred file
        :param file_base_name: base file name, e.g. Test_010101_HLA_A_Cov_bRP
        :param figname: particular name to append to the base name to speficy the figure.
        """
        image_dir = os.path.join(self._output_dir, "image_output")
        if not os.path.exists(image_dir):
            os.mkdir(image_dir)
        fig.savefig(os.path.join(image_dir, file_base_name + "_" + figname + ".png"))

    def init_logo_directory(self, file_base_name) -> str:
        """
        Initializes logo_results directory to store logo motif output from palmotif.
        :param file_base_name: the base name of the file under which to store logo results
        """
        logo_dir = os.path.join(self._output_dir, "logo_results")
        if not os.path.exists(logo_dir):
            os.mkdir(logo_dir)
        specific_logo_dir = os.path.join(logo_dir, file_base_name)  # Assumes each file has a unique name
        if not os.path.exists(specific_logo_dir):
            os.mkdir(specific_logo_dir)
        return specific_logo_dir

    def append_union_data(self, to_add: pd.DataFrame):
        """
        Appends data to the union table, aggregating duplicates by incrementing the Count column.
        :param to_add: the new dataframe to add
        """
        print("Union table current length:", len(self._union_table))
        self._union_table = pd.concat([self._union_table, to_add], ignore_index=True)
        print("Union table second length:", len(self._union_table))
        self._union_table = self._union_table.groupby(
            ['Allele', 'HLAP_sequence'], as_index=False, sort=False
        ).agg({
            'HLAP_master_descriptions': 'first',
            'HLAP_length': 'first',
            'Count': 'sum'
        })

    def append_overview_data(self, overview_data: dict):
        """
        Appends a new output row to the overview data
        """
        self._final_table = pd.concat([self._final_table, pd.DataFrame([overview_data])], ignore_index=True)

    def save_union_file(self):
        """
        Saves the current version of the union file
        """
        self._union_table.to_csv(self._union_table_path)

    def save_overview_file(self):
        """
        Saves the current version of the overview file
        """
        self._final_table.to_csv(self._final_table_path)
