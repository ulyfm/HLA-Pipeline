import os
import unittest

from hlapipeline.HLAPipeline import HLAPipeline
from hlapipeline.files.OutputAdapter import OutputAdapter
from hlapipeline.files.PeptideFileAdapter import PeptideFileAdapter


class HLAPipelineTests(unittest.TestCase):
    def setUp(self):
        self.peptide_file1 = PeptideFileAdapter(os.path.join("test_files", "test1_000000_HLA_A010101_HEL_bRP_PeptideGroups.txt"))
        output = OutputAdapter("test_output", os.path.join("test_output", "union_table.csv"),
                               os.path.join("test_output", "overview_table.csv"))
        self.pipeline1 = HLAPipeline(output, self.peptide_file1, None)

    '''
    File adapter tests
    '''

    def test_file_name(self):
        self.assertEqual(self.peptide_file1.get_file_name(), "test1_000000_HLA_A010101_HEL_bRP_PeptideGroups.txt")

    def test_base_name(self):
        self.assertEqual(self.peptide_file1.get_base_name(), "test1_000000_HLA_A010101_HEL_bRP")

    def test_date_created(self):
        self.assertEqual(self.peptide_file1.get_date_created(), "000000")

    def test_get_cotransduced(self):
        self.assertEqual(self.peptide_file1.get_cotransduced(), "HEL")

    def test_allele_inference(self):
        self.assertEqual(self.peptide_file1.get_allele(), 'HLA_A010101')
        self.assertEqual(self.pipeline1.final_row['HLA_allele'], 'HLA_A010101')

    def test_standardize_table(self):
        self.assertTrue("HLAP_sequence" in self.pipeline1.get_data_table().columns)
        self.assertTrue("HLAP_RT" in self.pipeline1.get_data_table().columns)
        self.assertTrue("HLAP_master_accessions" in self.pipeline1.get_data_table().columns)
        self.assertTrue("HLAP_master_descriptions" in self.pipeline1.get_data_table().columns)
        self.assertTrue("HLAP_length" in self.pipeline1.get_data_table().columns)

    def test_standardize_table_exceptions(self):
        self.assertEqual(True, False)  # add assertion here

    def test_standardize_infer_length(self):
        self.assertEqual(True, False)  # add assertion here

    '''
    HLA pipeline tests
    '''

    def test_remove_sp(self):
        self.pipeline1.remove_sp()
        self.assertTrue("AAAAGSLSR" not in self.pipeline1.get_data_table()['HLAP_sequence'])
        self.assertEqual(self.pipeline1.final_row['sp_count'], 1)

    def test_remove_frag(self):
        self.pipeline1.remove_frag()
        self.assertTrue("AAAAGGAFP" not in self.pipeline1.get_data_table()['HLAP_sequence'])
        self.assertEqual(self.pipeline1.final_row['fragment_count'], 1)

    def test_remove_duplicates(self):
        self.pipeline1.remove_duplicates()
        self.assertEqual(len(self.pipeline1.get_data_table()), len(self.pipeline1.get_data_table().drop_duplicates(subset="HLAP_sequence")))

    def test_identify_cotransduced(self):
        self.pipeline1.identify_cotransduced(peptide=".*iroquois.*")
        self.assertEqual(self.pipeline1.final_row['co-transduced_protein_0'], "Iroquois-class homeodomain protein IRX-4 OS=Homo sapiens OX=9606 GN=IRX4 PE=1 SV=2")
        self.assertEqual(self.pipeline1.final_row['co-transduced_peptides_0'], 'AAAAAATSLSQ')

    def test_identify_cotransduced_multiple(self):
        self.pipeline1.identify_cotransduced(peptide=".*LAM")
        self.assertEqual(self.pipeline1.final_row['co-transduced_protein_0'], "FLAM")
        self.assertEqual(self.pipeline1.final_row['co-transduced_peptides_0'], 'FLAMFLAM, GLAMGLAM')
        self.assertEqual(self.pipeline1.final_row['co-transduced_protein_1'], "BLAM")
        self.assertEqual(self.pipeline1.final_row['co-transduced_peptides_1'], 'CLAMFLAM')
        self.assertEqual(self.pipeline1.final_row['co-transduced_protein_2'], "ZLAM")
        self.assertEqual(self.pipeline1.final_row['co-transduced_peptides_2'], 'CLAMFLAM')
        self.assertEqual(self.pipeline1.final_row['co-transduced_protein_3'], "SLAM")
        self.assertEqual(self.pipeline1.final_row['co-transduced_peptides_3'], 'FLAMCLAM')
        self.assertEqual(self.pipeline1.final_row['co-transduced_protein_4'], "GLAM")
        self.assertEqual(self.pipeline1.final_row['co-transduced_peptides_4'], 'CLAMCLAM')

    def test_generate_logos(self):
        self.assertEqual(True, False)  # add assertion here

    def test_db_data(self):
        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
