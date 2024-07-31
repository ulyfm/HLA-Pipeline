import os
import numpy as np
import pandas as pd


def save_csv(table: pd.DataFrame, main_dir: str, folder_name: str, file_name: str) -> None:
    """
    Helper function to save a file according to the following format:

    Main output directory
    > Sub folder
    > > File

    :param table: the table to save
    :param main_dir: the main output dir
    :param folder_name: the sub folder to save in
    :param file_name: the file name to save as
    :return: None
    """
    if not os.path.exists(main_dir):
        os.mkdir(main_dir)
    if not folder_name == "":
        main_dir = os.path.join(main_dir, folder_name)
        if not os.path.exists(main_dir):
            os.mkdir(main_dir)
    table.to_csv(
        os.path.join(main_dir, file_name))


def get_flat_files(path_array: list, exclude_dirs: list, use_peptide_group_files: bool, base_path="") -> list:
    """
    Returns the appropriate files to process from the base directory, with basic filtering applied.

    :param path_array: a list of paths to search through.
    :param exclude_dirs: a list of directory names to exclude, for example if the output directory is the same as the input.
    :param use_peptide_group_files: if true, this means that only files that end in PeptideGroups.txt are included.
    :param base_path: the base path to (optionally) use a a relative directory that contains the paths in path_array
    """
    ret = []
    for path in path_array:
        apath = os.path.join(base_path, path)
        if not os.path.exists(apath):
            print("file not found:", apath)
            continue
        else:
            if path.startswith("."):
                print("ignored file:", path)
                continue
            if os.path.isdir(apath):
                if path in exclude_dirs:
                    print("ignored directory:", apath)
                    continue
                else:
                    ret += get_flat_files(os.listdir(apath), exclude_dirs, use_peptide_group_files, base_path=apath)
            else:
                if use_peptide_group_files:
                    if apath.endswith("PeptideGroups.txt"):
                        ret.append(apath)
                    else:
                        print("ignored file:", apath)
                else:
                    ret.append(apath)
    return ret
