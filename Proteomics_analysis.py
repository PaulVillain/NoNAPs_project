### Pre-processing of the data

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text


MS_data_jan = pd.read_csv('/Users/pvillain/Documents/Projects/NAPs_project/Nucleoid_proteomics/Proteomics_Oxford/January_2025_data/raw_data_jan25.csv')
MS_data_oct = pd.read_csv('/Users/pvillain/Documents/Projects/NAPs_project/Nucleoid_proteomics/Proteomics_Oxford/October_2024_data/MS_raw_data_oct.csv')
MS_data_jan = MS_data_jan.rename(columns={MS_data_jan.columns[0]: "protein_accessions_1"})
MS_data_oct = MS_data_oct.rename(columns={MS_data_oct.columns[0]: "protein_accessions_1"})
gff_IDs = pd.read_csv('/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/NextSeq/correspondence_Ecoli_strains/Manitoba_MG1655_corres_table_xy_no_duplicate.csv')
gff_IDs = gff_IDs.rename(columns={"protein_Manitoba": "protein_accessions_1"})

# Filter out duplicated rows based on the specified column
col_name = 'protein_accessions_1'
duplicated_rows = gff_IDs[gff_IDs.duplicated(subset=col_name, keep=False)]

# Get the list of names of duplicated values -> IS elements
duplicated_values = duplicated_rows[col_name].unique().tolist()

# Keep only the first occurrence of each duplicated value in the 'protein_accessions_1' column
gff_IDs = gff_IDs.drop_duplicates(subset='protein_accessions_1', keep='first')

# Merge the two datasets with the annotation table
MS_data_IDs_jan = gff_IDs.merge(MS_data_jan, how='inner', on='protein_accessions_1')
MS_data_IDs_jan.to_csv('/Users/pvillain/Documents/Projects/NAPs_project/Nucleoid_proteomics/Proteomics_Oxford/Analysis_2025/MS_data_IDs_jan.csv')
MS_data_IDs_oct = gff_IDs.merge(MS_data_oct, how='inner', on='protein_accessions_1')
MS_data_IDs_oct.to_csv('/Users/pvillain/Documents/Projects/NAPs_project/Nucleoid_proteomics/Proteomics_Oxford/Analysis_2025/MS_data_IDs_oct.csv')


# Merge the two datasets together
MS_data_IDs_Ox = MS_data_IDs_oct.merge(MS_data_IDs_jan, how='outer', on='protein_accessions_1')
MS_data_IDs_Ox.to_csv('/Users/pvillain/Documents/Projects/NAPs_project/Nucleoid_proteomics/Proteomics_Oxford/Analysis_2025/MS_data_IDs_Ox.csv')

# average the values by condition (IBAQ)

# WCE
MS_data_IDs_Ox['WT_exp_WCE_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_wce_wt_exp1', 'raw_iBAQ_wce_wt_exp3']].mean(axis=1)
MS_data_IDs_Ox['WT_stat_WCE_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_wce_wt_stat1', 'raw_iBAQ_wce_wt_stat3']].mean(axis=1)

MS_data_IDs_Ox['delta1_exp_WCE_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_wce_delta1_1', 'raw_iBAQ_wce_delta1_2']].mean(axis=1)

MS_data_IDs_Ox['delta2_exp_WCE_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_wce_delta2_1', 'raw_iBAQ_wce_delta2_2']].mean(axis=1)

MS_data_IDs_Ox['delta3_exp_WCE_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_wce_delta3_1', 'raw_iBAQ_wce_delta3_2']].mean(axis=1)

# 'raw_iBAQ_wce_delta4_1' replicate not considered -> contamination with delta1 strain
MS_data_IDs_Ox['delta4_exp_WCE_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_wce_delta4_2']].mean(axis=1)

MS_data_IDs_Ox['delta5_exp_WCE_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_wce_delta5_1', 'raw_iBAQ_wce_delta5_2']].mean(axis=1)

MS_data_IDs_Ox['delta6_exp_WCE_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_wce_delta6_1', 'raw_iBAQ_wce_delta6_2']].mean(axis=1)

MS_data_IDs_Ox['delta7_exp_WCE_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_wce_delta7_1', 'raw_iBAQ_wce_delta7_2']].mean(axis=1)

MS_data_IDs_Ox['delta8_exp_WCE_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_wce_delta8_1', 'raw_iBAQ_wce_delta8_2']].mean(axis=1)

MS_data_IDs_Ox['delta9_exp_WCE_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_wce_delta9_exp1', 'raw_iBAQ_wce_delta9_exp3']].mean(axis=1)
MS_data_IDs_Ox['delta9_stat_WCE_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_wce_delta9_stat1', 'raw_iBAQ_wce_delta9_stat3']].mean(axis=1)

# top fraction
MS_data_IDs_Ox['WT_exp_top_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_top_wt_exp1', 'raw_iBAQ_top_wt_exp3']].mean(axis=1)
MS_data_IDs_Ox['WT_stat_top_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_top_wt_stat1', 'raw_iBAQ_top_wt_stat3']].mean(axis=1)

# 'raw_iBAQ_top_delta4_1' replicate not considered -> contamination with delta1 strain
MS_data_IDs_Ox['delta4_exp_top_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_top_delta4_2']].mean(axis=1)

MS_data_IDs_Ox['delta5_exp_top_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_top_delta5_1', 'raw_iBAQ_top_delta5_2']].mean(axis=1)

MS_data_IDs_Ox['delta9_exp_top_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_top_delta9_exp1', 'raw_iBAQ_top_delta9_exp3']].mean(axis=1)
MS_data_IDs_Ox['delta9_stat_top_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_top_delta9_stat1', 'raw_iBAQ_top_delta9_stat3']].mean(axis=1)

# nucleoid fraction
MS_data_IDs_Ox['WT_exp_nucleoid_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_nuc_wt_exp1', 'raw_iBAQ_nuc_wt_exp3']].mean(axis=1)
MS_data_IDs_Ox['WT_stat_nucleoid_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_nuc_wt_stat1', 'raw_iBAQ_nuc_wt_stat3']].mean(axis=1)

# 'raw_iBAQ_nuc_delta4_1' replicate not considered -> contamination with delta1 strain
MS_data_IDs_Ox['delta4_exp_nucleoid_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_nuc_delta4_2']].mean(axis=1)

MS_data_IDs_Ox['delta5_exp_nucleoid_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_nuc_delta5_1', 'raw_iBAQ_nuc_delta5_2']].mean(axis=1)

MS_data_IDs_Ox['delta9_exp_nucleoid_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_nuc_delta9_exp1', 'raw_iBAQ_nuc_delta9_exp3']].mean(axis=1)
MS_data_IDs_Ox['delta9_stat_nucleoid_avg_IBAQ'] = MS_data_IDs_Ox[['raw_iBAQ_nuc_delta9_stat1', 'raw_iBAQ_nuc_delta9_stat3']].mean(axis=1)



# average the values by condition (raw intensities)

# WCE
MS_data_IDs_Ox['WT_exp_WCE_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_wce_wt_exp1', 'raw_quantity_wce_wt_exp3']].mean(axis=1)
MS_data_IDs_Ox['WT_stat_WCE_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_wce_wt_stat1', 'raw_quantity_wce_wt_stat3']].mean(axis=1)

MS_data_IDs_Ox['delta1_exp_WCE_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_wce_delta1_1', 'raw_quantity_wce_delta1_2']].mean(axis=1)

MS_data_IDs_Ox['delta2_exp_WCE_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_wce_delta2_1', 'raw_quantity_wce_delta2_2']].mean(axis=1)

MS_data_IDs_Ox['delta3_exp_WCE_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_wce_delta3_1', 'raw_quantity_wce_delta3_2']].mean(axis=1)

MS_data_IDs_Ox['delta4_exp_WCE_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_wce_delta4_1', 'raw_quantity_wce_delta4_2']].mean(axis=1)

MS_data_IDs_Ox['delta5_exp_WCE_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_wce_delta5_1', 'raw_quantity_wce_delta5_2']].mean(axis=1)

MS_data_IDs_Ox['delta6_exp_WCE_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_wce_delta6_1', 'raw_quantity_wce_delta6_2']].mean(axis=1)

MS_data_IDs_Ox['delta7_exp_WCE_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_wce_delta7_1', 'raw_quantity_wce_delta7_2']].mean(axis=1)

MS_data_IDs_Ox['delta8_exp_WCE_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_wce_delta8_1', 'raw_quantity_wce_delta8_2']].mean(axis=1)

MS_data_IDs_Ox['delta9_exp_WCE_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_wce_delta9_exp1', 'raw_quantity_wce_delta9_exp3']].mean(axis=1)
MS_data_IDs_Ox['delta9_stat_WCE_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_wce_delta9_stat1', 'raw_quantity_wce_delta9_stat3']].mean(axis=1)

# top fraction
MS_data_IDs_Ox['WT_exp_top_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_top_wt_exp1', 'raw_quantity_top_wt_exp3']].mean(axis=1)
MS_data_IDs_Ox['WT_stat_top_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_top_wt_stat1', 'raw_quantity_top_wt_stat3']].mean(axis=1)

MS_data_IDs_Ox['delta4_exp_top_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_top_delta4_1', 'raw_quantity_top_delta4_2']].mean(axis=1)

MS_data_IDs_Ox['delta5_exp_top_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_top_delta5_1', 'raw_quantity_top_delta5_2']].mean(axis=1)

MS_data_IDs_Ox['delta9_exp_top_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_top_delta9_exp1', 'raw_quantity_top_delta9_exp3']].mean(axis=1)
MS_data_IDs_Ox['delta9_stat_top_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_top_delta9_stat1', 'raw_quantity_top_delta9_stat3']].mean(axis=1)

# nucleoid fraction
MS_data_IDs_Ox['WT_exp_nucleoid_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_nuc_wt_exp1', 'raw_quantity_nuc_wt_exp3']].mean(axis=1)
MS_data_IDs_Ox['WT_stat_nucleoid_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_nuc_wt_stat1', 'raw_quantity_nuc_wt_stat3']].mean(axis=1)

MS_data_IDs_Ox['delta4_exp_nucleoid_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_nuc_delta4_1', 'raw_quantity_nuc_delta4_2']].mean(axis=1)

MS_data_IDs_Ox['delta5_exp_nucleoid_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_nuc_delta5_1', 'raw_quantity_nuc_delta5_2']].mean(axis=1)

MS_data_IDs_Ox['delta9_exp_nucleoid_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_nuc_delta9_exp1', 'raw_quantity_nuc_delta9_exp3']].mean(axis=1)
MS_data_IDs_Ox['delta9_stat_nucleoid_avg_raw_intensity'] = MS_data_IDs_Ox[['raw_quantity_nuc_delta9_stat1', 'raw_quantity_nuc_delta9_stat3']].mean(axis=1)


# Add the DBP column from proteinfer to tell if a protein is a DNA binding protein
df_ProteInfer = pd.read_csv("/Users/pvillain/Documents/Projects/NAPs_project/ProteInfer_BW25113/GCF_015534855_protein_orderedSup40_DBP.csv")
df_ProteInfer = df_ProteInfer.rename(columns={"sequence_name": "protein_accessions_1", "description": "DBP"})
df_ProteInfer = df_ProteInfer.drop(['predicted_label', 'confidence'], axis=1)
MS_data_IDs_Ox_DBP = pd.merge(MS_data_IDs_Ox, df_ProteInfer, on="protein_accessions_1", how='left')
# Replace empty strings with NaN and remove rows where all values are empty
MS_data_IDs_Ox_DBP.replace('', pd.NA, inplace=True)
MS_data_IDs_Ox_DBP = MS_data_IDs_Ox_DBP.dropna(how='all')
MS_data_IDs_Ox_DBP.to_csv('/Users/pvillain/Documents/Projects/NAPs_project/Nucleoid_proteomics/Proteomics_Oxford/Analysis_2025/MS_data_IDs_Ox_DBP.csv')



### Plot the protein abundance

# WT_exp_IBAQ values
# Load the DataFrame and preprocess
MS_data_IDs_DBP = pd.read_csv('/Users/pvillain/Documents/Projects/NAPs_project/Nucleoid_proteomics/Proteomics_Oxford/Analysis_2025/MS_data_IDs_Ox_DBP.csv')
MS_data_IDs_DBP['WT_exp_WCE_avg_IBAQ'] = pd.to_numeric(MS_data_IDs_DBP['WT_exp_WCE_avg_IBAQ'], errors='coerce')

# Fill NaN values with 0
MS_data_IDs_DBP.fillna(0, inplace=True)

# Sort the DataFrame by 'WT_exp_WCE_avg_IBAQ'
MS_data_IDs_DBP = MS_data_IDs_DBP.sort_values('WT_exp_WCE_avg_IBAQ')

# Add x_coord and y_coord columns to the existing DataFrame
MS_data_IDs_DBP['x_coord'] = range(1, len(MS_data_IDs_DBP) + 1)  # Rank (starting at 1)
MS_data_IDs_DBP['y_coord'] = MS_data_IDs_DBP['WT_exp_WCE_avg_IBAQ'].astype(int)  # Convert y to integer

# Scatter plot
ax = MS_data_IDs_DBP.plot.scatter(x='x_coord', y='y_coord', figsize=(8, 6), c='lightgray')

# List of NAPs
list_NAPs = ['hupA', 'hupB', 'hns', 'ihfB', 'ihfA', 'stpA', 'dps', 'fis', 'lrp']

# Annotate DNA binding proteins and adjust their positions
texts = []
for _, row in MS_data_IDs_DBP.iterrows():
    if row['DBP'] == 'DNA binding':
        color = 'purple'
        plt.scatter(row['x_coord'], row['y_coord'], c=color)
        if row['Gene Name'] in ['rpoZ', 'rpoA']:  # Combining conditions for rpoZ and rpoA
            texts.append(plt.text(row['x_coord'], row['y_coord'], row['Gene Name'], color='black', ha='center', va='bottom', fontsize=9))
for _, row in MS_data_IDs_DBP.iterrows():
    if row['Gene Name'] in list_NAPs:
        color = 'red'
        plt.scatter(row['x_coord'], row['y_coord'], c=color)
        texts.append(plt.text(row['x_coord'], row['y_coord'], row['Gene Name'], color='black', ha='right', va='bottom', fontsize=9))

# Scale and label adjustments
plt.yscale('log')
sns.despine(offset=15)
plt.xlabel("Abundance Rank")
plt.ylabel("log10 Copies per cell")
plt.title("NAPs over DNA binding proteins (WCE_WT_expo)")

# Adjust text positions manually
for text in texts:
    x, y = text.get_position()
    text.set_position((x + 10, y))  # Adjust the offset (e.g., 10)
adjust_text(texts)

# Save the plot
plt.savefig("/Users/pvillain/Documents/Projects/NAPs_project/Nucleoid_proteomics/Proteomics_Oxford/plots/2025/WCE_WT_expo_rank_abundance_with_other_NAPs_candidates.pdf", dpi=300, transparent=True)
plt.show()




### Plot nucleoid abundance

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, linregress
from adjustText import adjust_text

# Load data
MS = pd.read_csv(
    '/Users/pvillain/Documents/Projects/NAPs_project/'
    'Nucleoid_proteomics/Proteomics_Oxford/Analysis_2025/'
    'MS_data_IDs_Ox_DBP.csv'
)

# Columns we need
intensity_cols = [
    'WT_exp_nucleoid_avg_IBAQ',
    'delta9_exp_nucleoid_avg_IBAQ'
]


# ---- Plot ----
plt.figure(figsize=(8, 6))

# Scatter all points
plt.scatter(x, y, c='lightgray', alpha=0.2, s=15)

# x = y line

line_min = min(x.min(), y.min())
line_max = max(x.max(), y.max())
plt.plot(
    [line_min, line_max],
    [line_min, line_max],
    color='black',
    lw=1.5,
    label='x = y'
)

sns.despine(offset=15)
plt.xlabel("Protein abundance in WT nucleoid fraction (IBAQ)")
plt.ylabel("Protein abundance in ∆NAP9 nucleoid fraction (IBAQ)")
plt.title("Protein abundance in the nucleoid fraction")
plt.legend()
plt.xscale('log')
plt.yscale('log')
# plt.xlim(-0.1, 1E5)
# plt.ylim(-0.1, 1E5)

list_POI = ['gyrA']
texts = []

# DNA-binding proteins in purple
dbp_mask = MS.get('DBP', pd.Series(index=MS.index, dtype=object)).eq('DNA binding')
MS_dbp = MS.loc[dbp_mask]
plt.scatter(MS_dbp['WT_exp_nucleoid_avg_IBAQ'], MS_dbp['delta9_exp_nucleoid_avg_IBAQ'], c='purple', s=15)

# Protein of interest in red (and label them)
gene_col = 'Gene Name_x_y' if 'Gene Name_x_y' in MS.columns else 'Gene Name'
if gene_col in MS.columns:
    MS_POI = MS.loc[MS[gene_col].isin(list_POI)]
    if not MS_POI.empty:
        plt.scatter(MS_POI['WT_exp_nucleoid_avg_IBAQ'], MS_POI['delta9_exp_nucleoid_avg_IBAQ'], c='red', s=15)
        for _, r in MS_POI[['WT_exp_nucleoid_avg_IBAQ','delta9_exp_nucleoid_avg_IBAQ', gene_col]].iterrows():
            texts.append(
                plt.text(
                    r['WT_exp_nucleoid_avg_IBAQ'],
                    r['delta9_exp_nucleoid_avg_IBAQ'],
                    str(r[gene_col]),
                    color='black', fontsize=7,
                    ha='center', va='bottom'
                )
            )

# Adjust text positions
if texts:
    adjust_text(texts)

plt.tight_layout()
# plt.savefig(
#     '/Users/pvillain/Documents/Projects/NAPs_project/Nucleoid_proteomics/Proteomics_Oxford/plots/2025/'
#     'Nucleoid_abundance_WT_vs_NAP9_IBAQ_with_regression.pdf'
# )
plt.show()





### Plot nucleoid enrichment

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from adjustText import adjust_text

# Load data
MS = pd.read_csv(
    '/Users/pvillain/Documents/Projects/NAPs_project/'
    'Nucleoid_proteomics/Proteomics_Oxford/Analysis_2025/'
    'MS_data_IDs_Ox_DBP.csv'
)

# Columns we need
intensity_cols = [
    'WT_exp_nucleoid_avg_raw_intensity',
    'WT_exp_top_avg_raw_intensity',
    'delta9_exp_nucleoid_avg_raw_intensity',
    'delta9_exp_top_avg_raw_intensity'
]

# ---- Strict filtering pipeline ----
# 1) Remove existing inf -> NaN, then drop NaNs in required columns
MS = MS.replace([np.inf, -np.inf], np.nan)
before_n = len(MS)
MS = MS.dropna(subset=intensity_cols).copy()

# 2) Keep only rows where all required intensities are strictly > 0
MS = MS.loc[(MS[intensity_cols] > 0).all(axis=1)].copy()
after_intensity_filter = len(MS)
print(f"Removed {before_n - after_intensity_filter} rows due to null/inf or non-positive intensities; {after_intensity_filter} rows remain before ratio calc.")

# 3) Compute enrichment ratios
MS['enrichment_WT'] = MS['WT_exp_nucleoid_avg_IBAQ'] / MS['WT_exp_top_avg_IBAQ']
MS['enrichment_delta9'] = MS['delta9_exp_nucleoid_avg_IBAQ'] / MS['delta9_exp_top_avg_IBAQ']

# 4) Keep only rows with finite, positive ratios
ratio_mask = (
    np.isfinite(MS['enrichment_WT']) & np.isfinite(MS['enrichment_delta9']) &
    (MS['enrichment_WT'] > 0) & (MS['enrichment_delta9'] > 0)
)
MS = MS.loc[ratio_mask].copy()

# 5) Log2 transform
MS['log2_enrichment_WT'] = np.log2(MS['enrichment_WT'])
MS['log2_enrichment_delta9'] = np.log2(MS['enrichment_delta9'])

# 6) Final safety: finite logs only
finite_mask = np.isfinite(MS['log2_enrichment_WT']) & np.isfinite(MS['log2_enrichment_delta9'])
MS = MS.loc[finite_mask].copy()
print(f"Final dataset for plotting: {len(MS)} rows.")

# ---- Plot ----
plt.figure(figsize=(8, 6))

# Baseline scatter (all points in light gray)
plt.scatter(
    MS['log2_enrichment_WT'],
    MS['log2_enrichment_delta9'],
    c='lightgray', alpha=0.2, s=15
)

sns.despine(offset=15)
plt.xlabel("WT log₂ enrichment (nucleoid versus top fraction)")
plt.ylabel("NAP9 log₂ enrichment (nucleoid versus top fraction)")
plt.title("Nucleoid enrichment WT vs NAP9 (exponential phase)")

list_NAPs = ['crp', 'mutM']
texts = []

# 1) DNA-binding proteins in purple (vectorized)
dbp_mask = MS.get('DBP', pd.Series(index=MS.index, dtype=object)).eq('DNA binding')
MS_dbp = MS.loc[dbp_mask]
if not MS_dbp.empty:
    plt.scatter(
        MS_dbp['log2_enrichment_WT'],
        MS_dbp['log2_enrichment_delta9'],
        c='purple', s=15
    )

# 2) Known NAPs in red (and label them)
gene_col = 'Gene Name_x_y' if 'Gene Name_x_y' in MS.columns else 'Gene Name'
if gene_col in MS.columns:
    MS_naps = MS.loc[MS[gene_col].isin(list_NAPs)]
    if not MS_naps.empty:
        plt.scatter(
            MS_naps['log2_enrichment_WT'],
            MS_naps['log2_enrichment_delta9'],
            c='red', s=15
        )
        for _, r in MS_naps[['log2_enrichment_WT','log2_enrichment_delta9', gene_col]].iterrows():
            texts.append(
                plt.text(
                    r['log2_enrichment_WT'],
                    r['log2_enrichment_delta9'],
                    str(r[gene_col]),
                    color='black', fontsize=7,
                    ha='center', va='bottom'
                )
            )

# Nudge overlapping labels
if texts:
    adjust_text(texts)

plt.tight_layout()
# plt.savefig('/Users/pvillain/Documents/Projects/NAPs_project/Nucleoid_proteomics/Proteomics_Oxford/plots/2025/Nucleoid_enrichment_WT_vs_NAP9_(exponential_phase).pdf')
plt.show()






### Plot abundance of DNA-binding proteins (DBP) in whole cell extract (WCE)

import pandas as pd
import matplotlib.pyplot as plt

# Files
proteomics_file = "/Users/pvillain/Documents/Manuscripts/NAPs_paper/MS_data/processed/MS_data_IDs_DBP_prophages_2025.csv"
prophage_file = "/Users/pvillain/Documents/Projects/NAPs_project/Hi-C/reference_genomes/genes_in_prophages.csv"

# Load data
data = pd.read_csv(proteomics_file)
prophages = pd.read_csv(prophage_file)

# Columns
gene_col_proteomics = "Gene Name_x_y"
gene_col_prophage = "gene"

conditions = {
    'WT WCE expo': 'WT_exp_WCE_avg_IBAQ',
    'WT WCE stat': 'WT_stat_WCE_avg_IBAQ',
    'Δ9 WCE expo': 'delta9_exp_WCE_avg_IBAQ',
    'Δ9 WCE stat': 'delta9_stat_WCE_avg_IBAQ',
}

# Clean proteomics table
data = data.dropna(subset=[gene_col_proteomics])

for col in conditions.values():
    data[col] = pd.to_numeric(data[col], errors="coerce")

# Keep only protein-coding prophage genes
prophages = prophages[prophages["gene_biotype"] == "protein_coding"].copy()

# Define gene categories
rnap_genes = ['rpoA', 'rpoB', 'rpoC', 'rpoD', 'rpoZ']
topoisomerase_genes = ['gyrA', 'gyrB', 'topA', 'topB', 'parC', 'parE']
nap_genes = ['hupA', 'hupB', 'hns', 'ihfA', 'ihfB', 'fis', 'dps', 'lrp', 'stpA']

# ---- Category sums ----

category_sums = pd.DataFrame(
    index=[
        "Prophage DBP",
        "Other DNA-binding proteins",
        "RNA polymerase",
        "Topoisomerases",
        "NAPs"
    ],
    columns=list(conditions.values())
)

# Prophage genes
prophage_genes = set(prophages[gene_col_prophage])

# 1. Aggregate prophage genes that are also DNA-binding proteins
prophage_data = data[
    data[gene_col_proteomics].isin(prophage_genes) &
    (data["DBP"] == "DNA binding")
]

category_sums.loc["Prophage DBP", list(conditions.values())] = (
    prophage_data[list(conditions.values())].sum(skipna=True)
)

# 2. RNAP category
rnap_data = data[
    data[gene_col_proteomics].isin(rnap_genes)
]

category_sums.loc["RNA polymerase", list(conditions.values())] = (
    rnap_data[list(conditions.values())].sum(skipna=True)
)

# 3. Topoisomerases category
topoisomerase_data = data[
    data[gene_col_proteomics].isin(topoisomerase_genes)
]

category_sums.loc["Topoisomerases", list(conditions.values())] = (
    topoisomerase_data[list(conditions.values())].sum(skipna=True)
)

# 4. NAP category
nap_data = data[
    data[gene_col_proteomics].isin(nap_genes)
]

category_sums.loc["NAPs", list(conditions.values())] = (
    nap_data[list(conditions.values())].sum(skipna=True)
)

# 5. Other DNA-binding proteins
# Exclude RNAP genes, topoisomerases, NAPs, and all prophage genes so there is no double counting
dbp_data = data[
    (data["DBP"] == "DNA binding") &
    (~data[gene_col_proteomics].isin(rnap_genes)) &
    (~data[gene_col_proteomics].isin(topoisomerase_genes)) &
    (~data[gene_col_proteomics].isin(nap_genes)) &
    (~data[gene_col_proteomics].isin(prophage_genes))
]

category_sums.loc["Other DNA-binding proteins", list(conditions.values())] = (
    dbp_data[list(conditions.values())].sum(skipna=True)
)

# Convert to numeric
category_sums = category_sums.astype(float)

# Convert each category abundance to % of total proteome per condition
plot_table = pd.DataFrame(index=category_sums.index)

for label, col in conditions.items():
    total_proteome = data[col].sum(skipna=True)
    plot_table[label] = (category_sums[col] / total_proteome) * 100

# Optional: sort categories by total contribution
plot_table["total"] = plot_table.sum(axis=1)
plot_table = plot_table.sort_values("total", ascending=False).drop(columns="total")

# Plot
ax = plot_table.T.plot(
    kind="bar",
    stacked=True,
    figsize=(9, 8)
)

ax.set_ylabel("Proportion of total proteome (%)")
ax.set_title("Prophage DBP, other DNA-binding proteins, RNAP, topoisomerases and NAP abundance - IBAQ values")

plt.xticks(rotation=20, ha="right")
plt.legend(
    title="Category",
    bbox_to_anchor=(1.05, 1),
    loc="upper left"
)

plt.tight_layout()

# plt.savefig(
#     "/Users/pvillain/Documents/Projects/NAPs_project/Nucleoid_proteomics/Correlation_with_published_MS_dataset/stacked_barplot_DBP_percentage_full_internal_dataset.pdf",
#     format="pdf",
#     transparent=True,
#     bbox_inches="tight",
#     dpi=300
# )
plt.show()
