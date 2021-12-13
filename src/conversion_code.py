# Imports

from datetime import datetime
import os

from Bio import Entrez
import numpy as np
import pandas as pd
import pytz

import WoRMS  # custom functions for querying WoRMS API


def load_data(filename):
    """
    Loads a dataframe into pandas using read_csv()

    :param filename: file name (string)
    :return df: pandas data frame containing data from csv file
    """
    filename = os.getcwd().replace('src', os.path.join('raw', filename))
    df = pd.read_csv(filename)

    return df


def fill_lowest_taxon(df, col_names):
    """
    Takes the occurrence pandas data frame and fills missing values in scientificName
    with values from the first non-missing taxonomic rank column. The names of the taxonomic
    rank columns are listed in col_names.

    :param df: pandas data frame
    :param col_names: list of column names (list of strings)
    :return df: pandas data frame with scientificName populated
    """
    col_names.reverse()

    for col in col_names[:-1]:
        df['scientificName'] = df['scientificName'].combine_first(df[col])

    col_names.reverse()

    return df


def find_not_matched(df, name_dict):
    """
    Takes the occurrence pandas data frame and name_dict matching scientificName values
    with names on WoRMS and returns a list of names that did not match on WoRMS.

    :param df: pandas data frame
    :param name_dict: dictionary matching scientificName values with names matched on WoRMS
    :return unmatched: list of names that did not match on WoRMS (list of strings)
    """
    unmatched = []

    for n in df['scientificName'].unique():
        if n not in name_dict.keys():
            unmatched.append(n)

    try:
        unmatched.remove(np.nan)
    except ValueError:
        pass

    return unmatched


def replace_not_matched(df, unmatched, col_names):
    """
    Takes the occurrence pandas data frame and a list of scientificName values that
    did not match on WoRMS and replaces those values with NaN in the columns specified by col_names.

    :param df: pandas data frame
    :param unmatched: list of scientificName values that did not match on WoRMS (list of strings)
    :param col_names: list of column names (list of strings)
    :return df: pandas data frame with NaN instead of names that did not match on WoRMS
    """
    df[col_names] = df[col_names].replace(unmatched, np.nan)

    return df


# Load data
print('Loading data...')
plate = load_data('asv_table.csv')
meta = load_data('metadata_table.csv')

# Sequence_ID column in the plate dataframe uniquely identifies a water sample
print('Starting conversion...')
occ = pd.DataFrame({'eventID': plate['Sequence_ID']})

# Merge with meta to obtain columns that can be added directly from metadata
metadata_cols = [
    'seqID',
    'eventDate',
    'decimalLatitude',
    'decimalLongitude',
    'env_broad_scale',
    'env_local_scale',
    'env_medium',
    'target_gene',
    'primer_sequence_forward',
    'primer_sequence_reverse',
    'pcr_primer_name_forward',
    'pcr_primer_name_reverse',
    'pcr_primer_reference',
    'sop',
    'seq_meth',
    'samp_vol_we_dna_ext',
    'nucl_acid_ext',
    'nucl_acid_amp',
]

dwc_cols = metadata_cols.copy()
dwc_cols[0] = 'eventID'

occ = occ.merge(meta[metadata_cols], how='left', left_on='eventID', right_on='seqID')
occ.drop(columns='seqID', inplace=True)
occ.columns = dwc_cols

# Format eventDate
pst = pytz.timezone('America/Los_Angeles')
eventDate = [pst.localize(datetime.strptime(dt, '%m/%d/%y %H:%M')).isoformat() for dt in occ['eventDate']]
occ['eventDate'] = eventDate

# Clean columns
occ['sop'] = occ['sop'].str.replace('|', ' | ', regex=False)
occ = occ.rename(columns={'primer_sequence_forward': 'pcr_primer_forward',
                          'primer_sequence_reverse': 'pcr_primer_reverse'})

# Add extension terms that weren't in metadata file (obtained by asking data provider)
occ['target_subfragment'] = 'V9'
occ['lib_layout'] = 'paired'
occ['otu_class_appr'] = 'dada2;version;ASV'
occ['otu_seq_comp_appr'] = 'blast;version;80% identity | MEGAN6;version; bitscore:100:2%'  # Needs to be updated
occ['otu_db'] = 'Genbank nr;221'  # Needs to be updated

# Create an occurrenceID that will uniquely identify each ASV observed within a water sample
occ['occurrenceID'] = plate.groupby('Sequence_ID')['ASV'].cumcount()+1
occ['occurrenceID'] = occ['eventID'] + '_occ' + occ['occurrenceID'].astype(str)

# Add DNA_sequence
occ['DNA_sequence'] = plate['ASV']

# Add scientificName, taxonomic info
occ['scientificName'] = plate['Species']
occ['kingdom'] = plate['Kingdom']
occ['phylum'] = plate['Phylum']
occ['class'] = plate['Class']
occ['order'] = plate['Order']
occ['family'] = plate['Family']
occ['genus'] = plate['Genus']

# Replace 'unknown', 'unassigned', etc. in scientificName and taxonomy columns with NaN
cols = ['scientificName', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus']
occ[cols] = occ[cols].replace({'unassigned': np.nan,
                               's_': np.nan,
                               'g_': np.nan,
                               'unknown': np.nan,
                               'no_hit': np.nan})

# Get unique species names
names = occ['scientificName'].unique()
names = names[~pd.isnull(names)]  # remove NaN

# Replace non-Linnaean species names with NaN
non_latin_names = []
for name in names:
    words_in_name = name.split(' ')
    if len(words_in_name) > 2:
        non_latin_names.append(name)
non_latin_names_dict = {i: np.nan for i in non_latin_names}

non_latin_names_dict['phototrophic eukaryote'] = np.nan
non_latin_names_dict['Candida <clade Candida/Lodderomyces clade>'] = np.nan
occ['scientificName'].replace(non_latin_names_dict, inplace=True)

# Replace entries where kingdom = 'Eukaryota' with the WoRMS-approved 'Biota'
occ.loc[occ['kingdom'] == 'Eukaryota', 'kingdom'] = 'Biota'

# Iterate to match lowest possible taxonomic rank on WoRMS
# Note that cols (list of taxonomic column names) was defined in a previous code block
print('Attempting to match scientific names on WoRMS. This may take several minutes.')
name_name_dict = {}
name_id_dict = {}
name_taxid_dict = {}
name_class_dict = {}

not_matched = [1]

while len(not_matched) > 0:
    # Step 1 - fill
    occ = fill_lowest_taxon(occ, cols)

    # Step 2 - get names to match
    to_match = find_not_matched(occ, name_name_dict)

    # Step 3 - match
    print('Attempting to match {num} names on WoRMS...'.format(num=len(to_match)))
    name_id, name_name, name_taxid, name_class = WoRMS.run_get_worms_from_scientific_name(to_match, verbose_flag=False)
    name_id_dict = {**name_id_dict, **name_id}
    name_name_dict = {**name_name_dict, **name_name}
    name_taxid_dict = {**name_taxid_dict, **name_taxid}
    name_class_dict = {**name_class_dict, **name_class}

    # Step 4 - get names that didn't match
    not_matched = find_not_matched(occ, name_name_dict)
    print('Number of names not matched: {num}'.format(num=len(not_matched)))

    # Step 5 - replace these values with NaN
    occ = replace_not_matched(occ, not_matched, cols)

print('Matching complete.')

# Change scientificName to Biota in cases where all taxonomic information is missing
occ.loc[occ['scientificName'].isna() == True, 'scientificName'] = 'Biota'

# Fix taxonomy columns
occ[cols[1:]] = plate[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']].copy()
occ[cols[1:]] = occ[cols[1:]].replace({
    'unassigned': '',
    's_': '',
    'g_': '',
    'unknown': '',
    'no_hit': ''})

# Add scientific name-related columns
occ['scientificNameID'] = occ['scientificName']
occ['scientificNameID'].replace(name_id_dict, inplace=True)
occ['taxonID'] = occ['scientificName']
occ['taxonID'].replace(name_taxid_dict, inplace=True)
occ['scientificName'].replace(name_name_dict, inplace=True)
occ['nameAccordingTo'] = 'WoRMS'

# Get set up to query NCBI taxonomy
# ----- Insert your email here -----
Entrez.email = 'dianalg@mbari.org'
# ----------------------------------

# Get list of all databases available through this tool
record = Entrez.read(Entrez.einfo())
all_dbs = record['DbList']

# Get NCBI taxIDs for each name in dataset
print('Attempting to retrieve {} taxonomy IDs from NCBI...'.format(len(names)))
name_ncbiid_dict = {}
for name in names:
    handle = Entrez.esearch(db='taxonomy', retmax=10, term=name)
    record = Entrez.read(handle)
    name_ncbiid_dict[name] = record['IdList'][0]
    handle.close()

# Add NCBI taxonomy IDs under taxonConceptID
print('Completing conversion...')
name_ncbiid_dict['unassigned'], \
    name_ncbiid_dict['s_'], \
    name_ncbiid_dict['no_hit'], \
    name_ncbiid_dict['unknown'], \
    name_ncbiid_dict['g_'] = '', '', '', '', ''

occ['taxonConceptID'] = plate['Species'].copy()
occ['taxonConceptID'].replace(name_ncbiid_dict, inplace=True)

occ['taxonConceptID'] = 'NCBI:txid' + occ['taxonConceptID']
occ['taxonConceptID'].replace('NCBI:txid', '', inplace=True)

# Get identificationRemarks
occ = occ.merge(meta[['seqID', 'identificationRemarks']], how='left', left_on='eventID', right_on='seqID')
occ.drop(columns='seqID', inplace=True)

# Add name that matched in GenBank - i.e. the species name from the original data
occ['identificationRemarks'] = plate['Species'].copy() + ', ' + occ['identificationRemarks']

# Add basisOfRecord
occ['basisOfRecord'] = 'MaterialSample'

# Add identificationReferences
occ = occ.merge(meta[['seqID', 'identificationReferences']], how='left', left_on='eventID', right_on='seqID')
occ.drop(columns='seqID', inplace=True)
occ['identificationReferences'] = occ['identificationReferences'].str.replace('| ', ' | ', regex=False)

# Add organismQuantity (number of reads)
occ['organismQuantity'] = plate['Reads']
occ['organismQuantityType'] = 'DNA sequence reads'

# Add sampleSizeValue
count_by_seq = plate.groupby('Sequence_ID', as_index=False)['Reads'].sum()
occ = occ.merge(count_by_seq, how='left', left_on='eventID', right_on='Sequence_ID')
occ.drop(columns='Sequence_ID', inplace=True)
occ.rename(columns={'Reads': 'sampleSizeValue'}, inplace=True)

occ['sampleSizeUnit'] = 'DNA sequence reads'

# Add associatedSequences
occ = occ.merge(meta[['seqID', 'associatedSequences']], how='left', left_on='eventID', right_on='seqID')
occ.drop(columns='seqID', inplace=True)

# Drop records where organismQuantity = 0 (absences are not meaningful for this data set)
occ = occ[occ['organismQuantity'] > 0]

# Divide into occurrence and DNADerivedDataExt
ddd_cols = [
    'eventID',
    'occurrenceID',
    'DNA_sequence',
    'sop',
    'nucl_acid_ext',
    'samp_vol_we_dna_ext',
    'nucl_acid_amp',
    'target_gene',
    'target_subfragment',
    'lib_layout',
    'pcr_primer_forward',
    'pcr_primer_reverse',
    'pcr_primer_name_forward',
    'pcr_primer_name_reverse',
    'pcr_primer_reference',
    'seq_meth',
    'otu_class_appr',
    'otu_seq_comp_appr',
    'otu_db',
    'env_broad_scale',
    'env_local_scale',
    'env_medium',
]

DNADerivedData = occ[ddd_cols].copy()
occ.drop(ddd_cols[2:], axis=1, inplace=True)

# Save
print('Saving data...')

folder = os.getcwd().replace('src', 'processed')
occ_filename = os.path.join(folder, 'occurrence.csv')
ddd_filename = os.path.join(folder, 'dna_extension.csv')

if not os.path.exists(folder):
    os.makedirs(folder)

occ.to_csv(occ_filename, index=False, na_rep='NaN')
DNADerivedData.to_csv(ddd_filename, index=False, na_rep='NaN')

print('Conversion complete.')
