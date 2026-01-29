import argparse
import sys

import pandas as pd
import os

from hlapipeline import LocalDatabase
from hlapipeline import HLAPipeline
from hlapipeline import util
from hlapipeline.files.OutputAdapter import OutputAdapter
from hlapipeline.files.PeptideFileAdapter import PeptideFileAdapter


def main():
    """
     Command line arguments. Descriptions for the arguments are in readme.md.
    """
    # Start by parsing all the command line arguments:
    parser = argparse.ArgumentParser(prog='HLA-Pipeline',
                                     description='Data pipeline for HLA mass spec processing')
    parser.add_argument('-o', '--output', default="hla_files")
    parser.add_argument('-i', '--input', default="hla_files")
    parser.add_argument('-u', '--union', default="union_table.csv")
    parser.add_argument('-t', '--overview', default="final_table.csv")
    parser.add_argument('-bulkcsv')
    parser.add_argument('-bulkdir')
    parser.add_argument('--peptide', action=argparse.BooleanOptionalAction)
    parser.add_argument('--skipcleanup', action=argparse.BooleanOptionalAction)
    parser.add_argument('--skipcotransduced', action=argparse.BooleanOptionalAction)
    parser.add_argument('--skipdataviz', action=argparse.BooleanOptionalAction)
    parser.add_argument('-d', '--database', action=argparse.BooleanOptionalAction)
    parser.add_argument('--updatedb', action=argparse.BooleanOptionalAction)
    parser.add_argument('--dbloc', action=argparse.BooleanOptionalAction)
    parser.add_argument('--assumecotransduced', action=argparse.BooleanOptionalAction)
    args = parser.parse_args()

    # Download (or just load) 3rd party databases.
    local_db = None
    if args.database or args.dbloc:
        local_db = LocalDatabase.LocalDatabase(args.updatedb, args.dbloc)

    # Next, find the data files we are going to be processing (although some may be excluded later by error checking)
    specific_files = [args.input]  # Directory containing files to access.
    exclude_dirs = ['sp_peptides', 'spRM_peptides', 'frag_peptides', 'fragRM_peptides', 'dup_peptides',
                    'final_peptides', 'final_peptides_8-14',
                    'image_output', 'logo_results', 'dupRM_fragRM_spRM_peptides']
    data_files_flat = util.get_flat_files(specific_files, exclude_dirs, not args.peptide)

    # If a final table exists in the output directory, we are going to append rows to it rather than creating a new one.
    outputdir = args.output
    unionfile = args.union
    overviewfile = args.overview
    if unionfile == 'union_table.csv':
        unionfile = os.path.join(outputdir, "union_table.csv")
    if overviewfile == 'final_table.csv':
        overviewfile = os.path.join(outputdir, "final_table.csv")
    output = OutputAdapter(outputdir, unionfile, overviewfile)

    if args.bulkcsv is None or args.bulkdir is None:
        """
            Main processing of HLA files:
        """
        for path in data_files_flat:
            print("Processing", path)
            peptide_file = PeptideFileAdapter(path)
            pipeline = HLAPipeline.HLAPipeline(output, peptide_file, local_db)
            print("File loaded.")
            # Generate final table output (i.e. metadata).
            result = pipeline.get_result(args.skipcleanup, args.skipcotransduced, args.assumecotransduced,
                                         args.skipdataviz, database_features=args.database)
            output.append_overview_data(result)
            # Generate complete output for inclusion in union table
            detailed_table = pipeline.get_data_table()[['HLAP_sequence', 'HLAP_length', 'HLAP_master_accessions',
                                                        'HLAP_master_descriptions']]
            detailed_table['Count'] = 1
            print("HLA allele:", result['HLA_allele'])
            detailed_table['Allele'] = result['HLA_allele']
            # Merge results with union table via aggregation. (Count increases with multiple copies)
            output.append_union_data(detailed_table)
            print("File finished processing.")
            print()
        output.save_union_file()
        print("Union file saved")
        output.save_overview_file()
        print("Overview file stored as final_table.csv")
    else:
        """
            Bulk processing of HLA files:
        """
        print("Bulk processing mode")
        bulk_dir = args.bulkdir
        bulk_csv_path = args.bulkcsv
        if not os.path.isdir(bulk_dir):
            print("Invalid bulk directory:", bulk_dir)
            sys.exit()
        if not os.path.exists(bulk_csv_path):
            print("Invalid bulk csv")
            sys.exit()
        bulk_csv = pd.read_csv(args.bulkcsv)
        overall_output = OutputAdapter(args.output, args.union, args.overview)
        for path in os.listdir(bulk_dir):
            if not os.path.isdir(os.path.join(bulk_dir, path)) or path in exclude_dirs:
                continue
            subset_output = OutputAdapter(os.path.join(args.output, path),
                                          os.path.join(args.output, path, path + "_union.csv"),
                                          os.path.join(args.output, path, path + "_overview.csv"))
            print("Processing allele directory:", path)
            for filepath in util.get_flat_files([os.path.join(bulk_dir, path)],
                                                exclude_dirs, not args.peptide):
                print("Processing", filepath)
                cotrans_peptides = list(bulk_csv[bulk_csv['file_name'] == filepath]['co-transduced protein(s)'])[0]
                print("with cotransduced:", cotrans_peptides)
                pipeline = HLAPipeline.HLAPipeline(subset_output, PeptideFileAdapter(path), local_db)
                print("File loaded.")
                # Generate output.
                result = pipeline.get_result(args.skipcleanup, args.skipcotransduced, False,
                                             args.skipdataviz, cotransduced_peptide=cotrans_peptides)
                subset_output.append_overview_data(result)
                overall_output.append_overview_data(result)
                detailed_table = pipeline.get_data_table()[['sequence', 'length', 'Master Protein Accessions',
                                                            'Master Protein Descriptions']]
                detailed_table['Count'] = 1
                detailed_table['Allele'] = result['HLA_allele']
                # Merge results with union table via aggregation. (Count increases with multiple copies)
                subset_output.append_union_data(detailed_table)
                overall_output.append_union_data(detailed_table)
                print("File finished processing.")
                print()
            subset_output.save_union_file()
            subset_output.save_overview_file()
            print(path, "allele union table saved.")
        overall_output.save_union_file()
        print("Overall union table saved.")
        overall_output.save_overview_file()
        print("Data stored as final_table.csv")


if __name__ == '__main__':
    main()
