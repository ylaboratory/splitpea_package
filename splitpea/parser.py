import pandas as pd
import re
from collections import defaultdict
import os 
import tempfile
import numpy as np  
import warnings 
import logging 

def process_suppa2(psivec_path, dpsi_path, map_path, splicing_events_filter=None, species="human", save=False, verbose=False):
    """
    Process and merge two splicing-related files (PSI and dPSI/p-val) into a unified DataFrame.

    Parameters:
      psivec_path (str): File path for the PSI file.
      dpsi_path (str): File path for the dPSI (and p-val) file.
      splicing_events_filter (list, optional): List of splicing event types to keep (e.g. ["SE", "A3"]).
                                               If None, no filtering is applied.
      species (str): Species name for gene symbol mapping (default is "human").

    Returns:
      final_df (pd.DataFrame): DataFrame with metadata columns and computed columns:
         - ensembl.id, symbol, chr, strand, exon.start, exon.end
         - psi columns (averaged per event group)
         - delta.psi and pval (from the dPSI file)
    """
    
    # Load the files using the first column as the index (row names)
    psivec_df = pd.read_csv(psivec_path, sep="\t", header=0, index_col=0)
    dpsi_df   = pd.read_csv(dpsi_path, sep="\t", header=0, index_col=0)
    
    # Merge the DataFrames on their index
    merged_df = psivec_df.merge(dpsi_df, left_index=True, right_index=True, 
                                suffixes=('_psivec', '_dpsi'))
    
    # Helper function to parse a row name into metadata.
    def parse_index(idx):
        """
        Parses a row name string assumed to be in the format:
          "ENSEMBL_ID;EVENT_INFO"
        where EVENT_INFO is colon-delimited, for example:
          "SE:X:76886222-76886503:76886613-76888625:-"
        Returns a Series with:
          ensembl.id, event, chr, strand, exon.start, exon.end.
        """
        try:
            parts = idx.split(";")
            ensembl_id = parts[0]
            event_info = parts[1]
            fields = event_info.split(":")
            event_type = fields[0]      
            chr_val    = fields[1]
            exon_range = fields[2].split("-")
            exon_start = exon_range[1] 
            exon_range2 = fields[3].split("-")
            exon_end   = exon_range2[0]
            strand     = fields[4]
        except Exception:
            ensembl_id, event_type, chr_val, strand, exon_start, exon_end = (np.nan,)*6
        return pd.Series([ensembl_id, event_type, chr_val, strand, exon_start, exon_end],
                         index=["ensembl.id", "event", "chr", "strand", "exon.start", "exon.end"])
    
    # Reset the index to parse the row names, then restore it
    merged_df = merged_df.reset_index()
    merged_df[['ensembl.id', 'event', 'chr', 'strand', 'exon.start', 'exon.end']] = merged_df['index'].apply(parse_index)
    merged_df = merged_df.set_index('index')
    
    # Optionally filter by splicing event types
    if splicing_events_filter is not None:
        merged_df = merged_df[ merged_df['event'].isin(splicing_events_filter) ]
    
    # # Retrieve gene symbols using mygene
    # mg = mygene.MyGeneInfo()
    # unique_ids = merged_df['ensembl.id'].dropna().unique().tolist()
    # results = mg.querymany(unique_ids, scopes='ensembl.gene', fields='symbol',
    #                        species=species, as_dataframe=True, verbose = verbose)
    # # Build mapping, ensuring that missing mappings result in np.nan.
    # mapping = results['symbol'].to_dict()
    # merged_df['symbol'] = merged_df['ensembl.id'].apply(lambda x: mapping.get(x, np.nan))

    map_df = pd.read_csv(map_path, sep='\t', usecols=['symbol', 'ensembl'])
    mapping = dict(zip(map_df['ensembl'], map_df['symbol']))
    merged_df['symbol'] = merged_df['ensembl.id'].map(mapping)
    
    all_cols = merged_df.columns.tolist()
    meta_cols_ = {
        'ensembl.id',
        'event',
        'chr',
        'strand',
        'exon.start',
        'exon.end',
        'symbol'
    }

    event_cols = [
        col for col in all_cols
        if "dpsi" not in col.lower()
        and "p-val" not in col.lower()
        and "pval" not in col.lower()
        and col not in meta_cols_
    ]
    
    # Clean column names by removing any file-specific suffixes (_psivec, _dpsi)
    def clean_col(col):
        return col.replace("_psivec", "").replace("_dpsi", "")
    event_cols_clean = [clean_col(col) for col in event_cols]
    
    # Group event columns by their base name by stripping any trailing underscore and digits.
    def base_group_name(col):
        return re.sub(r'_\d+$', '', col)
    
    groups = defaultdict(list)
    for orig_col, clean_name in zip(event_cols, event_cols_clean):
        base_name = base_group_name(clean_name)
        groups[base_name].append(orig_col)
    
    # Compute the mean for each event group.
    psi_cols = {}
    if len(groups) == 2:
        sorted_keys = sorted(groups.keys())
        psi_cols["psi." + str(sorted_keys[0])] = merged_df[groups[sorted_keys[0]]].mean(axis=1)
        psi_cols["psi." + str(sorted_keys[1])] = merged_df[groups[sorted_keys[1]]].mean(axis=1)
    else:
        print("WARNING THERE ARE MORE THAN 2 CONDITIONS")
        for group, cols in groups.items():
            psi_cols[f"psi.{group}"] = merged_df[cols].mean(axis=1)
    
    # Add the computed psi columns to the merged dataframe
    for new_col, series in psi_cols.items():
        merged_df[new_col] = series
    
    # Identify the dPSI and p-val columns by searching (ignoring case)
    dpsi_candidates = [col for col in merged_df.columns if "dpsi" in col.lower()]
    pval_candidates = [col for col in merged_df.columns if ("p-val" in col.lower() or "pval" in col.lower())]
    
    if dpsi_candidates:
        merged_df = merged_df.rename(columns={dpsi_candidates[0]: "delta.psi"})
    if pval_candidates:
        merged_df = merged_df.rename(columns={pval_candidates[0]: "pval"})
    
    # Define the final set of columns:
    # metadata, psi columns, then delta.psi and pval if available.
    meta_cols = ['ensembl.id', 'symbol', 'chr', 'strand', 'exon.start', 'exon.end']
    psi_final_cols = [col for col in merged_df.columns if col.startswith("psi.")]
    other_cols = []
    if "delta.psi" in merged_df.columns:
        other_cols.append("delta.psi")
    if "pval" in merged_df.columns:
        other_cols.append("pval")
    
    final_cols = meta_cols + psi_final_cols + other_cols
    final_df = merged_df[final_cols]
    final_df.loc[:, 'strand'] = final_df['strand'].map({'+': '1', '-': '-1'})
    final_df.loc[:, 'chr'] = final_df['chr'].str.replace(r'^chr', '', regex=True)


    #print(final_df.head(20))
    return final_df


def parse_suppa2(psivec_path, dpsi_path, map_path, splicing_events_filter=["SE"], species="human", verbose=False):
    """
    Parse SUPPA2 files using process_suppa2 and write the resulting DataFrame to a temporary file.
    The temporary file is written without row names, and it is saved in the same directory as the dpsi file.
    
    Returns:
      temp_file_path (str): Path to the temporary file containing the processed data.
    """
    final_df = process_suppa2(psivec_path, dpsi_path, map_path, splicing_events_filter, species, verbose=verbose)
    final_df['strand'] = final_df['strand'].replace({'+': '1', '-': '-1'})
    
    temp_dir = os.path.dirname(os.path.abspath(dpsi_path))
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', newline='', suffix=".txt", dir=temp_dir)
    final_df.to_csv(temp_file.name, index=False, sep="\t", na_rep="nan")
    temp_file.close()
    return temp_file.name



def calculate_psi_conds(events):
    '''
    Calculates psi for conditions 1 and 2 from rMATS file(s). Returns the mean of the events
    for each condition, or nan if all events are NA.

    Input:
    - events: a list with two inner lists, the first containing condition 1 values, the second
    containing condition 2 values
    '''
    events[0] = [np.nan if x == 'NA' else float(x) for x in events[0]]
    events[1] = [np.nan if x == 'NA' else float(x) for x in events[1]]

    psi_values_0 = [x for x in events[0] if str(x) != 'nan']
    psi_values_1 = [x for x in events[1] if str(x) != 'nan']

    dpsi_val_0 = -1
    dpsi_val_1 = -1

    if len(psi_values_0) == 0:
        dpsi_val_0 = 'nan'
    
    if len(psi_values_1) == 0:
        dpsi_val_1 = 'nan'
    
    # code sourced from SUPPA github
    # Ignore empty slice warning when calculating the mean
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', r'Mean of empty slice')
        if dpsi_val_0 != 'nan':
            dpsi_val_0 = np.nanmean(psi_values_0)

        if dpsi_val_1 != 'nan':
            dpsi_val_1 = np.nanmean(psi_values_1)
        
    return dpsi_val_0, dpsi_val_1 # mean, or nan if all values are NA


def parse_rmats(rmats_filepath, verbose=False): 
    '''
    Parses files from rMATS output. Generates output file in same directory as rmats_filepath.

    Input:
    - rmats_filepath: path to rMATS SE output, either JC.txt or JCEC.txt
    '''
    # place output file in same directory as .txt file
    rmats_dir = rmats_filepath.split('/')
    output_dir = ""
    for peice in range(len(rmats_dir)-1):
        output_dir += rmats_dir[peice] + "/"

    temp_dir = os.path.dirname(os.path.abspath(rmats_filepath))
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', newline='', suffix=".txt", dir=temp_dir)
    output = temp_file
    output.write('ensembl.id\tsymbol\tchr\tstrand\texon.start\texon.end\tpsi.condition1\tpsi.condition2\tdelta.psi\tpval\n')

    with open(rmats_filepath, 'r') as rmats_file:
        rmats_lines = rmats_file.readlines()

        for i in range(1, len(rmats_lines)): 
            rmats_contents = rmats_lines[i].strip().split('\t')
                    
            gene_id = rmats_contents[1].strip('"')  # Removes leading and trailing quotes
            gene_symbol = rmats_contents[2].strip('"')  # Removes leading and trailing quotes

            
            chr = rmats_contents[3].replace("chr", "", 1) # remove chr prefix
                    
            strand = "1"
            if rmats_contents[4] == '-': # - or + strand
                strand = "-1"
            
            start = rmats_contents[5] # exonStart_0base
            end = rmats_contents[6] #exonEnd

            cond1_values = rmats_contents[-3].split(",") # avg IncLevel1
            cond2_values = rmats_contents[-2].split(",") # avg IncLevel2
            events = [cond1_values, cond2_values]
            
            cond1, cond2 = calculate_psi_conds(events)

            if cond1 == 'nan' or cond2 == 'nan':
                dpsi = 'nan'
            else:
                dpsi = rmats_contents[-1] # IncLevelDifference

            pval = rmats_contents[-5] # PValue
            
            # filter out events where dpsi is nan
            if dpsi != 'nan':
                s = gene_id + "\t" + gene_symbol + "\t" + chr + "\t" + strand + "\t" + start + "\t" + end + "\t" + str(cond1) + "\t" + str(cond2) + "\t" + dpsi + "\t" + pval + "\n"
                output.write(s)
                    
    output.close()
    return temp_file.name