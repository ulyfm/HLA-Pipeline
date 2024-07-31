# HLA Pipeline

<!-- TOC -->
* [HLA Pipeline](#hla-pipeline)
  * [Overview](#overview)
  * [Installation](#installation)
    * [First-time Install](#first-time-install)
      * [Apple Silicon](#apple-silicon)
    * [Updating](#updating)
  * [Usage](#usage)
    * [1. Loading relevant files](#1-loading-relevant-files)
    * [2. Additional filtering](#2-additional-filtering)
    * [3. Co-transduced identification](#3-co-transduced-identification)
    * [4. Graph output](#4-graph-output)
    * [5. Append database data](#5-append-database-data)
    * [6. Overview and union table output](#6-overview-and-union-table-output)
    * [Bulk mode](#bulk-mode-)
<!-- TOC -->

## Overview

HLA Pipeline is a small program to assist with [ARTEMIS](https://doi.org/10.3389/fimmu.2021.658372), a platform for identifying
peptides bound by certain HLAs under varying cell conditions, such as co-transduction with a disease protein.

The functions that HLA Pipeline provides are:
- **Filtering:** Additional basic filtering of [ProteomeDiscoverer](https://www.thermofisher.com/us/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/proteome-discoverer-software.html) output data
- **Charts:** Automatic generation of peptide logo charts and length graphs
- **Tables:** Overview and union data tables summarizing output of multiple ARTEMIS runs
- **Identify Co-transduced:** Co-transduced peptide data mining
- (Optional) **Append Databases:** pulling additional protein data from several databases

## Installation

### First-time Install
Start by downloading the program. To do this, open Terminal and navigate to the folder where you want it to be installed.
Then clone the repository:
```commandline
git clone https://github.com/ulyfm/HLA-Pipeline.git
```
You may need to install command line tools if you are on mac (an automatic install popup should give you the option when you try to run git)

On windows, I strongly recommend [git bash](https://gitforwindows.org) as a way of running all of these commands. Another option is WSL.
The program is designed to be cross-compatible.

With the files downloaded, move inside the new folder and run the installation script. You can do that by typing these commands into Terminal:
```commandline
cd HLA-Pipeline
pip3 install .
```
Now you should be able to use the HLA-Pipeline tool from any directory by writing `hla-pipeline`.

#### Apple Silicon

There is an installation problem on new Apple computers. The installation of some of the dependencies
may get stuck on:
```
Building wheels for collected packages: parasail
Building wheel for parasail (setup.py)
```
The workaround is this:
1. Install [Homebrew](https://brew.sh)
2. Install automake: `brew install automake`
3. Install libtool: `brew install libtool`
4. Try to run the package installation again (`pip3 install .` from the HLA-Pipeline directory)

### Updating

*If you already have cloned the repository and want the most recent updates*, from inside the HLA-Pipeline directory execute the command:
```commandline
git pull
pip3 install .
```
To update the locally stored databases, use the `--updatedb` argument the next time you run the program.

## Usage

As discussed earlier, the program takes PeptideGroups output files from ProteomeDiscoverer and further processes them.

The basic process is:
1. Loading the relevant files
2. Additional filtering (duplicates, fragments, and 'sp' peptides)
3. Co-transduced identification
4. Graph output
5. (Optional) appending 3rd party database data
6. Overview and union table output

### 1. Loading relevant files
To start, you need to have a ProteomeDiscoverer PeptideGroups.txt output file. In a typical run,
you will put that file in a directory, then run `hla-pipeline` from the command line while it's your
working directory. For example, you might have a folder called `hla_files` on your desktop.

You would run:
```commandline
cd ~/Desktop
hla-pipeline
```

You may want to customize the working directory, or use separate ones for input and output. In that case, you can use
the following two command-line arguments:

**Set input directory/file:** (default: hla_files)
```commandline
hla-pipeline -i input_directory/file
```
**Set output directory:** (default: hla_files)
```commandline
hla-pipeline -o output_directory
```

So a typical run might look like:
```commandline
cd ~/Desktop
hla-pipeline -i hla_input -o hla_output
```

By default, the program only processes files ending in _PeptideGroups.txt, since this is what we expect from ProteomeDiscoverer output. In fact, any tab-separated file will work.
You can allow files with other endings by using the following command-line argument:

**Set whether to exclude files that don't end in _PeptideGroups.txt:** (default: True, so adding this argument makes it False)
```commandline
hla-pipeline --peptide
```
n.b. you must be careful with which files are in your input directory if you choose this option, since all of them
(except hidden files) will be included.

### 2. Additional filtering

'sp' (a.k.a. crap-ome) peptides, fragment peptides, and duplicates, are filtered by default in the respective order they were just listed.
At each stage,  additional files are saved for posterity. RM means removed, so spRM means 'sp' peptides have been removed
whereas 'sp_peptides' are all the peptides labeled 'sp' in their Master Accessions column.
- 'sp' (junk peptide) removal: spRM_peptides and sp_peptides.
- frag removal: fragRM_peptides and frag_peptides.
- duplicate removal: dupRM_fragRM_spRM_peptides and dup_peptides

n.b. removal is cumulative, so fragRM_peptides are a subset of spRM_peptides and dupRM_fragRM_spRM_peptides are a subset of
fragRM_peptides. (Sorry about the naming.)

This is generally good to include, but you can turn it off with the following command-line argument:

**Skip cleanup steps** (default: False.) This skips the cleanup steps (defrag, duplicate removal, junk peptide removal)
```commandline
hla-pipeline --skipcleanup
```

### 3. Co-transduced identification

At this stage, the program lets users identify co-transduced peptides as one of the processing steps.
After being prompted like so:
```
Co-transduced peptides - comma separated? (enter to skip)
>
```
You may enter a regular expression like `.*HLA.*`. The `.*` means any number of characters. So this expression matches
strings that contain HLA anywhere.

Note a couple of things:
- This is case-insensitive
- Both Master Protein Accessions and Master Protein Descriptions columns are compared against the regular expression
- You can have multiple peptides to be matched, simply by separating them with a comma like`.*HLA_A2.*, .*Histone`
- The characters ` `, `_`, and `-` are ignored by default unless they are part of your regular expression. For example,
`.*HLAA2.*` will match all strings that contain `HLA_A2`, `HLAA2`, `HLA A2`, etc. But if you specifically write out
the underscore, then `.*HLA_A2.*` will only match strings that contain `HLA_A2`.

To test your regular expressions, I recommend using [RegexPal](https://www.regexpal.com/).

The most useful symbol is going to be `.*`, which means any number (0+) of wildcard characters.

Some basic examples are here:
- Matching all strings beginning with "HLA": `^HLA.*`. `^` means it has to be at the beginning.
- Matching all strings including "HLA": `.*HLA.*`
- Matching all strings ending with "HLA": `.*HLA$`. `$` means it has to be at the ending.

n.b. only full matches are accepted.

There are two relevant command-line arguments for this stage:

**Skip co-transduced** (default: False, i.e. co-transduced input will be requested normally). This skips the step where
co-transduced peptides are queried.
```commandline
hla-pipeline --skipcotransduced
```

If you want to run a lot of co-transduced proteins in a row, Bulk Mode is recommended.

**Assume co-transduced** (default: False). This option makes the program assume co-transduced peptides from the file name.
A maximum of 1 peptide will be assumed and searched for. It will still work even with the 'skip co-transduced' option, although
that option will be overriden by assuming co-transduced (for the 1st protein only).
```commandline
hla-pipeline --assumecotransduced
```

The required format for assuming co-transduced is:
```
name_date_HLA_allele_cotransduced_protein_name_bRP_PeptideGroups.txt
```

The format depends on the first few underscores and where 'bRP' is, so you can put anything that doesn't contain
underscores for the other fields. It just needs to start with the fifth string and end with 'bRP'.

### 4. Graph output
This step outputs motif logo graphs and pie charts of peptide length.

You can skip it with this argument:

**Skip data visualization** (default: False, i.e. data viz will be generated normally). This skips the step where data
viz (pie charts and logo graphs) are generated.
```commandline
hla-pipeline --skipdataviz
```

### 5. Append database data

There is built-in functionality to append data from [Gene Ontology](https://geneontology.org)
and [ProteinAtlas](https://www.proteinatlas.org). This Includes the following data:

Gene Ontology:
- Gene name ('GO Name')
- Protein name ('GO Protein Name')
- Cellular components ('GO Cellular Components')

ProteinAtlas:
- Accession ('Accession_PA')
- Gene ('Gene')
- 'Gene description'
- 'Subcellular location'
- 'RNA tissue specificity'
- 'RNA tissue distribution'
- 'Evidence'
- 'RNA tissue specific nTPM'
- 'Tissue expression cluster'
- nTPM in HEK293 cells: 'HEK293 nTPM'
- 'TAU score - Tissue'
- 'Biological process'
- 'Molecular function'
- 'RNA single cell type specific nTPM'
- 'RNA tissue cell type enrichment'

Along with several other fields used to correlate proteins and genes across different databases.
The program must also download a Uniprot database that correlates protein identifiers across different
categories and databases.

To opt in to using the database features, use the following command:

**Use database features** (default: False). This opts in to appending extra 3rd party database data.
```commandline
hla-pipeline -d or --database
```

The program downloads the 3rd party databases for maximum speed (although it might take a while for the initial download).
Periodically, you may want to run the
following command to update the databases.

**Update local database** (default: False). This force updates the local database from the internet.
```commandline
hla-pipeline --updatedb
```

If you want to manually edit the databases being used, you can use this command to find where they are being stored
on your computer (it's OS specific):

**Database location** (default: False). This makes the database location print in console, then ends the program.
```commandline
hla-pipeline --dbloc
```

'go_cc_table.csv' is the summary table. To re-download a specific database, delete the database file as well as
'go_cc_table.csv' and it will be regenerated.

### 6. Overview and union table output

An 'overview' file (final_table.csv) with key information on co-transduced peptides, the number of peptides of
each length, and other details of the run is output.

A union table with key details of the peptides is also output. If a previous union file with the same name
exists, then they will be combined and the 'Count' column of overlapping peptides will be incremented.

You can customize the paths of these with the following two commands:

**Set union table:** (default: [output dir]/union_table.csv)
```commandline
hla-pipeline -u folder/union_table_path.csv
```
**Set overview table:** (default: [output dir]/final_table.csv)
```commandline
hla-pipeline -t folder/overview_table_path.csv
```

### Bulk mode ###

This is a new and exciting feature. Bulk mode is used like follows:
```commandline
hla-pipeline -bulkcsv bulk_csv -bulkdir bulk_file_directory
```

The files are expected to be in a shared directory (the "bulk file directory") separated into sub-directories
for each group of files you want a union table to be generated for. For example:

```
bulk_dir
 / allele1
   / file1.txt
   / file2.txt
 / allele2
   / file3.txt
```

Union tables will be placed in bulk_dir, as well as one for each allele.

The other component is the "bulk csv," which provides information for processing the files. The format is
as follows, with three columns ('HLA_allele','file_name','co-transduced protein(s)'). 'HLA_allele' is
currently not used. The program references 'file_name' to get the appropriate co-transduced peptides.
Here is an example, corresponding to the file structure above:

```
HLA_allele,file_name,co-transduced protein(s)
allele1,file1.txt,peptide1
allele1,file2.txt,peptide2
allele2,file3.txt,none
```