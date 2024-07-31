import os
import urllib.request
import pandas as pd
import platformdirs


class LocalDatabase:
    appname = "HLA-Pipeline"

    def __init__(self, force_update: bool, db_loc: bool):
        """
        Initializes local db access (table stored in memory for convenience/speed). If a cached file exists, no
        downloads will occur. Otherwise, the updated files will automatically be downloaded and parsed.
        :param force_update: if true, files will be replaced even if they exist. (this should be done every so often...)
        :param db_loc: whether to print the location of the database files AND THEN EXIT THE ENTIRE PROGRAM (!)
        """
        self._force_update = force_update
        user_data_folder = platformdirs.user_data_dir(self.appname)
        if not os.path.exists(user_data_folder):
            os.mkdir(user_data_folder)
        self._db_folder = os.path.join(user_data_folder, 'go_db_local')
        if not os.path.exists(self._db_folder):
            os.mkdir(self._db_folder)
        if db_loc:
            print(os.path.abspath(os.path.join(self._db_folder, "go_cc_table.csv")))
            exit(0)
        if not os.path.exists(os.path.join(self._db_folder, "go_cc_table.csv")) or force_update:
            self._go_table = pd.DataFrame()
            self.update_db()
        else:
            self._go_table = pd.read_csv(os.path.join(self._db_folder, "go_cc_table.csv"), index_col=0)

    def update_db(self):
        """
        Updates the database with files from the internet. The new data is cached on disk and used by other methods
        in this class until updated again.
        """
        if not os.path.exists(self._db_folder):
            os.mkdir(self._db_folder)
        # Start with every Uniprot accession code.
        self._go_table = self._parse_uniprot_idmapping()
        # Combine Gene Ontology terms and annotations.
        goterms = self._parse_go_obo()
        go_annotations = self._parse_go_associations()
        go_annotations = go_annotations.merge(goterms, on='GO_ID')
        go_annotations["GO Cellular Components"] = go_annotations["Qualifier"] + ":" + go_annotations["GO Term"]
        go_annotations = go_annotations.drop_duplicates()
        go_annotations = go_annotations.groupby(['Accession', 'GO Name', 'GO Protein Name'], as_index=False)[
            'GO Cellular Components'].agg(lambda x: "; ".join(x))

        self._go_table = self._go_table.merge(go_annotations, left_on="Accession", right_on="Accession", how="left")
        pa_table = self._parse_protein_atlas()
        pa_table = pa_table.astype(str)

        ensembl_subset = self._go_table[self._go_table['Ensembl'].notna()]
        ensembl_subset = ensembl_subset.merge(pa_table, on='Ensembl', how="left",
                                              suffixes=("", "_PA"))  # just removing duplicates...
        no_ensembl_subset = self._go_table[self._go_table['Ensembl'].isna()]
        no_ensembl_subset = no_ensembl_subset.merge(pa_table, on='Accession', how="left",
                                                    suffixes=("", "_PA"))  # just removing duplicates...
        self._go_table = pd.concat([ensembl_subset, no_ensembl_subset], ignore_index=True)
        print(self._go_table[self._go_table['Accession'] == 'L0R8F8'].to_string())
        pa_rna_table = self._parse_pa_hek293_ntpm()
        # Design choice: for nTPM, some protein atlas results can have multiple (why?) so we pick the max...
        pa_rna_table = pa_rna_table.groupby(['Gene name'], as_index=False).agg({
            'nTPM': 'max'
        })
        self._go_table = self._go_table.merge(pa_rna_table, left_on='Gene', right_on='Gene name', how="left")
        self._go_table = self._go_table.drop('Gene name', axis=1)

        self._go_table = self._go_table.rename({"nTPM": "HEK293 nTPM"}, axis=1)
        tau_data = self._parse_tau()
        self._go_table = self._go_table.merge(tau_data, on="Gene", how="left")
        print("Local databases finished loading.")

        def _join_no_na(x):
            ret = []
            for el in x:
                if not pd.isna(el):
                    ret.append(el)
            return "; ".join(ret)

        # To eliminate duplicate accessions, all data is merged semicolon-separated.
        self._go_table = self._go_table.groupby(['Accession'], as_index=False).agg({
            'Ensembl': _join_no_na,
            'GO Name': _join_no_na,
            'GO Protein Name': _join_no_na,
            'GO Cellular Components': _join_no_na,
            'Accession_PA': _join_no_na,
            'Gene': _join_no_na,
            'Gene description': _join_no_na,
            'Subcellular location': _join_no_na,
            'RNA tissue specificity': _join_no_na,
            'RNA tissue distribution': _join_no_na,
            'Evidence': _join_no_na,
            'RNA tissue specific nTPM': _join_no_na,
            'Tissue expression cluster': _join_no_na,
            'HEK293 nTPM': _join_no_na,
            'TAU score - Tissue': _join_no_na,
            'Biological process': _join_no_na,
            'Molecular function': _join_no_na,
            'RNA single cell type specific nTPM': _join_no_na,
            'RNA tissue cell type enrichment': _join_no_na
        })

        self._go_table.to_csv(os.path.join(self._db_folder, "go_cc_table.csv"))
        print("Local database: local files updated.")

    def _parse_go_obo(self) -> pd.DataFrame:
        """
        Download and parse gene ontology obo file, which connects gene ontology ids to phrases/attributes.
        :return: parsed data as a data frame.
        """
        print("Local database: downloading go.obo")
        if self._force_update or not os.path.exists(os.path.join(self._db_folder, "go.obo")):
            urllib.request.urlretrieve("http://purl.obolibrary.org/obo/go.obo",
                                       os.path.join(self._db_folder, "go.obo"))

        goterms = pd.DataFrame(columns=['GO_ID', 'GO Term'])
        index = 0
        print("Local database: parsing go.obo")
        with open(os.path.join(self._db_folder, "go.obo"), "rt") as obofile:
            term = False
            term_id = ''
            term_name = ''
            for line in obofile:
                if line.rstrip() == '[Term]':
                    if term:
                        term_id = ''
                        term_name = ''
                    term = True
                elif term and line.startswith("id: "):
                    term_id = line.rstrip()[4:]
                elif term and line.startswith("name: "):
                    term_name = line.rstrip()[6:]
                elif term and line.startswith("namespace: "):
                    term_namespace = line.rstrip()[11:]
                    if term_namespace != 'cellular_component':
                        term = False
                        term_id = ''
                        term_name = ''
                    else:
                        goterms.loc[index] = [term_id, term_name]
                        index += 1
        # os.remove(os.path.join(self._db_folder, "go.obo"))
        return goterms

    def _parse_go_associations(self) -> pd.DataFrame:
        """
        Download and parse gene ontology associations data ("cellular component" data only)
        For more info: https://geneontology.org/docs/ontology-documentation/
        :return: parsed data as a data frame.
        """
        print("Local database: downloading goa_human.gaf.gz")
        if self._force_update or not os.path.exists(os.path.join(self._db_folder, "goa_human.gaf.gz")):
            urllib.request.urlretrieve("http://geneontology.org/gene-associations/goa_human.gaf.gz",
                                       os.path.join(self._db_folder, "goa_human.gaf.gz"))
        print("Local database: parsing goa_human.gaf.gz")
        go_annotations = pd.read_csv(os.path.join(self._db_folder, "goa_human.gaf.gz"),
                                     compression='gzip',
                                     sep='\t',
                                     comment='!',
                                     dtype=str,
                                     names=['DB', 'Accession', 'GO Name', 'Qualifier', 'GO_ID',
                                            'DB:Reference', 'Evidence Code', 'With (or) From', 'Aspect',
                                            'GO Protein Name', 'DB Object Synonym', 'DB Object Type', 'Taxon',
                                            'Date', 'Assigned By', 'Annotation Extension', 'Gene Product Form ID'])
        go_annotations = go_annotations[['Accession', 'GO Name', 'Qualifier', 'GO_ID', 'GO Protein Name']]
        go_annotations = go_annotations[(go_annotations['Qualifier'] == 'located_in')
                                        | (go_annotations['Qualifier'] == 'part_of')
                                        | (go_annotations['Qualifier'] == 'is_active_in')
                                        | (go_annotations['Qualifier'] == 'colocalizes_with')]
        # os.remove(os.path.join(self._db_folder, "goa_human.gaf.gz"))
        return go_annotations

    def _parse_protein_atlas(self) -> pd.DataFrame:
        """
        Download and parse protein atlas data (subcellular location, tissue specificity, etc.)
        :return: parsed data as a data frame
        """
        print("Local database: downloading proteinatlas.tsv.zip")
        if self._force_update or not os.path.exists(os.path.join(self._db_folder, "proteinatlas.tsv.zip")):
            urllib.request.urlretrieve("https://www.proteinatlas.org/download/proteinatlas.tsv.zip",
                                       os.path.join(self._db_folder, "proteinatlas.tsv.zip"))
        print("Local database: parsing proteinatlas.tsv.zip")
        pa_table = pd.read_csv(os.path.join(self._db_folder, "proteinatlas.tsv.zip"),
                               compression='zip',
                               sep='\t',
                               dtype=str)
        pa_table = pa_table[['Uniprot', 'Gene', 'Gene description', 'Ensembl', 'Subcellular location',
                             'RNA tissue specificity', 'RNA tissue distribution', 'Evidence',
                             'RNA tissue specific nTPM', 'Tissue expression cluster', 'Biological process',
                             'Molecular function', 'RNA single cell type specific nTPM',
                             'RNA tissue cell type enrichment']]
        print("PA table columns:", pa_table.columns)
        pa_table = pa_table.rename(columns={'Uniprot': 'Accession'})
        # os.remove(os.path.join(self._db_folder, "proteinatlas.tsv.zip"))
        return pa_table

    def _parse_pa_hek293_ntpm(self) -> pd.DataFrame:
        """
        Download and parse ProteinAtlas HEK293 RNASeq gene expression data (units: nTPM)
        :return parsed data as a date frame
        """
        print("Local database: downloading rna_celline.tsv.zip")
        if self._force_update or not os.path.exists(os.path.join(self._db_folder, "rna_celline.tsv.zip")):
            urllib.request.urlretrieve("https://www.proteinatlas.org/download/rna_celline.tsv.zip",
                                       os.path.join(self._db_folder, "rna_celline.tsv.zip"))
        print("Local database: parsing rna_celline.tsv.zip")
        pa_rna_table = pd.read_csv(os.path.join(self._db_folder, "rna_celline.tsv.zip"),
                                   compression='zip',
                                   sep='\t',
                                   dtype=str)
        pa_rna_table = pa_rna_table[['Gene name', 'nTPM', 'Cell line']]
        pa_rna_table = pa_rna_table[pa_rna_table['Cell line'] == 'HEK293']
        # os.remove(os.path.join(self._db_folder, "rna_celline.tsv.zip"))
        return pa_rna_table

    def _parse_uniprot_idmapping(self) -> pd.DataFrame:
        """
        Download and parse UniProt ID mapping file. This is helpful because the relation of proteins to genes
        is n:1. (But some databases like ProteinAtlas only include 1 accession per gene...)
        :return parsed data as a date frame
        """
        print("Local database: downloading HUMAN_9606_idmapping_selected.tab.gz")
        if self._force_update or not os.path.exists(
                os.path.join(self._db_folder, "HUMAN_9606_idmapping_selected.tab.gz")):
            urllib.request.urlretrieve("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase"
                                       "/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz",
                                       os.path.join(self._db_folder, "HUMAN_9606_idmapping_selected.tab.gz"))
        print("Local database: parsing HUMAN_9606_idmapping_selected.tab.gz")
        idmapping = pd.read_csv(
            os.path.join(self._db_folder, "HUMAN_9606_idmapping_selected.tab.gz"),
            compression='gzip',
            sep='\t',
            dtype=str,
            header=None,
            usecols=[0, 18]
        )
        idmapping = idmapping.rename(columns={0: 'Accession', 18: 'Ensembl'})

        def _split_clean(ensembl: str):
            if pd.isna(ensembl):
                return pd.NA
            # print(ensembl)
            ar = ensembl.split(';')
            for i in range(0, len(ar)):
                if '.' in ar[i]:
                    ar[i] = ar[i][0:ar[i].index(".")].strip()
                else:
                    ar[i] = ar[i].strip()
            return ar

        idmapping['Ensembl'] = idmapping['Ensembl'].apply(_split_clean)
        idmapping = idmapping.explode('Ensembl').reset_index(drop=True)
        return idmapping

    def _parse_tau(self) -> pd.DataFrame:
        """
        Download and parse ProteinAtlas Tau data. Tau is a measure of tissue specificity of a particular gene.
        Tau of 1 means it is highly specific to a tissue type, whereas 0 means it shows up in all tissue types.
        :return parsed data as a date frame
        """
        print("Local database: downloading taudata.tsv")
        if self._force_update or not os.path.exists(os.path.join(self._db_folder, "taudata.tsv")):
            urllib.request.urlretrieve(
                "https://www.proteinatlas.org/api/search_download.php?search=&columns=g,"
                "t_RNA__tau&compress=no&format=tsv",
                os.path.join(self._db_folder, "taudata.tsv"))
        print("Local database: parsing taudata.tsv")
        tau_data = pd.read_csv(
            os.path.join(self._db_folder, "taudata.tsv"),
            sep='\t',
            dtype=str
        )
        tau_data = tau_data.drop_duplicates(subset="Gene") # TODO: is this necessary? what should the tiebreaker be?
        return tau_data

    def get_table(self):
        """
        :return: Gene ontology/protein atlas data table with relevant columns.
        """
        return self._go_table
