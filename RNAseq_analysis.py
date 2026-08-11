### Trimming and QC (TrimGalore)

#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=8:mem=50gb
#PBS -e trimm_error.txt
#PBS -o trimm_output.txt

module load anaconda3/personal

source activate /rds/general/user/pvillain/home/anaconda3/envs/RNAseq-env

cd /rds/general/user/pvillain/home/RNAseq/NAPs_project/231205_VH00504_171_AACMJ7KM5/Unaligned

for i in {1..28}; do
    trim_galore --cores 4 --fastqc --paired PV${i}_R1.fastq PV${i}_R2.fastq
done

trim_galore --cores 4 --fastqc --paired Undetermined_R1.fastq Undetermined_R2.fastq


qsub HPC_trim_galore_loop.sh



### rRNA in silico depletion (Bowtie2)

# align the reads on rRNA and retrieve both reads that maps to rRNA (aligned_rRNA) and reads that do not match (rRNA_clean)

#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=30:mem=50gb
#PBS -e bowtie2_RNAseq_error.txt
#PBS -o bowtie2_RNAseq_output.txt

module load anaconda3/personal

source activate /rds/general/user/pvillain/home/anaconda3/envs/bowtie2

/rds/general/user/pvillain/home/RNAseq/NAPs_project

for i in {1..28}; do
    bowtie2 --threads 30 \
        -x /rds/general/user/pvillain/home/RNAseq/NAPs_project/genomes/rRNA_BW25113 \
        -1 /rds/general/user/pvillain/home/RNAseq/NAPs_project/231205_VH00504_171_AACMJ7KM5/Unaligned/fq_after_trimming/PV$i\_R1_val_1.fq \
        -2 /rds/general/user/pvillain/home/RNAseq/NAPs_project/231205_VH00504_171_AACMJ7KM5/Unaligned/fq_after_trimming/PV$i\_R2_val_2.fq \
        -S /rds/general/user/pvillain/home/RNAseq/NAPs_project/rRNA_aligned/PV$i\_aligned_rRNA.sam \
        --al-conc /rds/general/user/pvillain/home/RNAseq/NAPs_project/rRNA_aligned/PV$i\_aligned_rRNA.fq \
        --un-conc /rds/general/user/pvillain/home/RNAseq/NAPs_project/rRNA_clean/PV$i\_rRNAclean.fq; \
done

qsub HPC_bowtie2_RNAseq_loop_rRNA.sh



### Align the reads on E. coli and B. subtilis genomes (spike-in) (Bowtie2)

#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=30:mem=50gb
#PBS -e bowtie2_RNAseq_error_2.txt
#PBS -o bowtie2_RNAseq_output_2.txt

module load anaconda3/personal

source activate /rds/general/user/pvillain/home/anaconda3/envs/bowtie2

cd /rds/general/user/pvillain/home/RNAseq/NAPs_project

for i in {1..28}; do bowtie2 -x /rds/general/user/pvillain/home/RNAseq/NAPs_project/genomes/Bsubtilis_168 -1 /rds/general/user/pvillain/home/RNAseq/NAPs_project/rRNA_clean/PV$i\_rRNAclean.1.fq -2 /rds/general/user/pvillain/home/RNAseq/NAPs_project/rRNA_clean/PV$i\_rRNAclean.2.fq -S /rds/general/user/pvillain/home/RNAseq/NAPs_project/Bsub_aligned/PV$i\_aligned_Bsub.sam; done
for i in {1..28}; do bowtie2 -x /rds/general/user/pvillain/home/RNAseq/NAPs_project/genomes/BW25113 -1 /rds/general/user/pvillain/home/RNAseq/NAPs_project/rRNA_clean/PV$i\_rRNAclean.1.fq -2 /rds/general/user/pvillain/home/RNAseq/NAPs_project/rRNA_clean/PV$i\_rRNAclean.2.fq -S /rds/general/user/pvillain/home/RNAseq/NAPs_project/Ecoli_aligned/PV$i\_aligned_Ecoli.sam; done


qsub HPC_bowtie2_RNAseq_loop_Bsub_Ecoli.sh 



### Index and sort the sam files (SamTools)

#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=10:mem=50gb
#PBS -e samtools_reads_error.txt


module load anaconda3/personal

source activate /rds/general/user/pvillain/home/anaconda3/envs/samtools-env

cd /rds/general/user/pvillain/home/RNAseq/NAPs_project/Ecoli_aligned

for i in {1..28}; do samtools view -@ 6 -c -F 260 /rds/general/user/pvillain/home/RNAseq/NAPs_project/Ecoli_aligned/PV$i\_aligned_Ecoli.sam >> counts_Ecoli.txt; done
for i in {1..28}; do samtools view -@ 6 -c -F 260 /rds/general/user/pvillain/home/RNAseq/NAPs_project/Bsub_aligned/PV$i\_aligned_Bsub.sam >> counts_Bsub.txt; done

for i in {1..28}; do samtools sort -@ 10 /rds/general/user/pvillain/home/RNAseq/NAPs_project/Ecoli_aligned/PV$i\_aligned_Ecoli.sam -o /rds/general/user/pvillain/home/RNAseq/NAPs_project/Ecoli_aligned/PV$i\_aligned_Ecoli_srt.bam; done
for i in {1..28}; do samtools sort -@ 10 /rds/general/user/pvillain/home/RNAseq/NAPs_project/Bsub_aligned/PV$i\_aligned_Bsub.sam -o /rds/general/user/pvillain/home/RNAseq/NAPs_project/Bsub_aligned/PV$i\_aligned_Bsub_srt.bam; done

for i in {1..28}; do samtools index -@ 10 /rds/general/user/pvillain/home/RNAseq/NAPs_project/Ecoli_aligned/PV$i\_aligned_Ecoli_srt.bam; done
for i in {1..28}; do samtools index -@ 10 /rds/general/user/pvillain/home/RNAseq/NAPs_project/Bsub_aligned/PV$i\_aligned_Bsub_srt.bam; done

qsub HPC_samtools_count_sort_index.sh




### Make the counts on E. coli and B. subtilis genomes (HTseq)

#!/bin/bash
for i in {1..28}
do
   qsub -v index="${i}" HPC_HTseq_NAPs_RS.sh
done

# script HPC_HTseq_NAPs_RS.sh:

#!/bin/bash
#PBS -l walltime=01:30:00
#PBS -l select=1:ncpus=1:mem=5gb

module load anaconda3/personal

source activate /rds/general/user/pvillain/home/anaconda3/envs/HTseq-env

cd /rds/general/user/pvillain/home/RNAseq/NAPs_project

htseq-count -n 1 -t gene -r pos -s no -f bam \
	--idattr=locus_tag \
	--additional-attr=locus_tag \
	--additional-attr=gene \
	--additional-attr=Name \
	--nonunique all \
	--stranded=no \
	"/rds/general/user/pvillain/home/RNAseq/NAPs_project/Ecoli_aligned/Ecoli_srt_bam/PV${index}_aligned_Ecoli_srt.bam" \
	"/rds/general/user/pvillain/home/RNAseq/NAPs_project/genomes/GCF_015534855.1_genomic.gff" \
	> "/rds/general/user/pvillain/home/RNAseq/NAPs_project/Ecoli_aligned/Ecoli_HTseq_count_2/PV${index}_Ecoli_HTseq.tsv"

htseq-count -n 1 -t gene -r pos -s no -f bam \
	--idattr=locus_tag \
	--additional-attr=locus_tag \
	--additional-attr=gene \
	--additional-attr=Name \
	--nonunique all \
	--stranded=no \
	"/rds/general/user/pvillain/home/RNAseq/NAPs_project/Bsub_aligned/Bsub_srt_bam/PV${index}_aligned_Bsub_srt.bam" \
	"/rds/general/user/pvillain/home/RNAseq/NAPs_project/genomes/GCF_000009045.1_ASM904v1_genomic.gff" \
	> "/rds/general/user/pvillain/home/RNAseq/NAPs_project/Bsub_aligned/Bsub_HTseq_count_2/PV${index}_Bsub_HTseq.tsv"


qsub HPC_HTseq_parallel.sh




### Differential analysis (DEseq2)

#Formatting of the DEseq2 input

import os
import pandas as pd

# Specify the directory containing the '.tsv' files
dossier = "/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/NextSeq/HTseq_output"

# Iterate through all files in the directory
for nom_fichier in os.listdir(dossier):
    # Check if the file has the '.tsv' extension
    if nom_fichier.endswith('.tsv'):
        # Create the full path to the file
        chemin_fichier = os.path.join(dossier, nom_fichier)
        
        # Read the content of the '.tsv' file into a DataFrame
        df = pd.read_csv(chemin_fichier, sep='\t', header=None, names=["gene", "description", "name", "count"])
        
        # Drop the last five rows from the DataFrame
        df = df.iloc[:-5]
        
        # Save the updated DataFrame back to the same file with a '.txt' extension
        new_fichier = os.path.splitext(chemin_fichier)[0] + '.txt'
        df.to_csv(new_fichier, sep='\t', index=False)

import pandas as pd 
import glob 
import os

#create an empty dataframe
df = pd.DataFrame(columns=['gene'])

#loop through all the TXT files in the directory
for file in glob.glob("*.txt"):
    # get the file name without the extension
    file_name = os.path.splitext(file)[0]
    # read each TXT file as a pandas dataframe
    file_df = pd.read_csv(file, sep="\t")
    # rename the 'count' column with the file name
    file_df.rename(columns={'count': file_name}, inplace=True)
    # merge the dataframes by the 'gene' column
    df = pd.merge(df, file_df[['gene', file_name]], on='gene', how='outer')

df.set_index('gene', inplace=True)

#fill missing values with 0
df.fillna(0, inplace=True)

#save the merged dataframe as a CSV file
df.to_csv('merged_data.csv')

df = pd.read_csv('merged_data.csv')

#transpose the dataframe to fit with DEseq2 requirements
df = df.transpose()

#Create the second input, the "column" dataframe
#directly with excel -> merged_data_2.csv and 
#get rid of the header and index
df2 = pd.read_csv('merged_sampled_data_2.csv', index_col=0)


# DEseq2 analysis (RStudio)


# Install DEseq2
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  install.packages("htmltools")
  # Dowload and install gfortran package
  BiocManager::install("DESeq2")

# Libraries
library(DESeq2)
  
# Load metaData
metaData <- read.csv('/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/NextSeq/HTseq_output/coldata.csv', header = TRUE, sep = ",")
# Change metaData$condition (remove last number)
metaData$condition=sapply(metaData$condition, function(x) gsub(pattern = "\\d+$", replacement = "", x))
# Attribute rownames
rownames(metaData)=metaData$X
metaData=metaData[ , -1]
# Change class of metaData$condition to be factors
metaData$condition=as.factor(metaData$condition)

# Load countData
countData <- read.csv('/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/NextSeq/HTseq_output/merged_data.csv', header = TRUE, sep = ",")
head(countData)
# Change column names in countData
New_names=sapply(colnames(countData), function(x) gsub(pattern = "_.*", replacement = "", x))
colnames(countData)=New_names
# Sort countData columns to be in the same order than metaData rownames
New_order_columns=c("gene", rownames(metaData))
countData_sorted=countData[ , New_order_columns]

# Run DEseq2
dds <- DESeqDataSetFromMatrix(countData=countData_sorted, 
                              colData=metaData, 
                              design=~condition, tidy = TRUE)

dds2 <- DESeq(dds)

# Take a look at the results table
res <- results(dds2,contrast = c("condition","9E","9S"),pAdjustMethod = "fdr")
head(res)

head(results(dds2, tidy=TRUE)) #let's look at the results table


# Summary of differential gene expression
summary(res) #summary of results




### Replicates comparison by PCA

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# import the data as pandas dataframe
os.chdir("/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/NextSeq/figures")

df = pd.read_csv(
    "/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/Precise_1K/comparison_Precise_NAPs/figures_input/PCA_logTPM_NAPs_samples.csv",
    index_col=0
)

# Transpose so that rows = samples, columns = genes
df_t = df.T

# Standardize features (genes) across samples
scaler = StandardScaler()
df_st = scaler.fit_transform(df_t)

# Perform PCA and get sample coordinates
pca = PCA()
scores = pca.fit_transform(df_st)

# Explained variance
explained_variance_ratio = pca.explained_variance_ratio_

print("Explained Variance for each Principal Component:")
for i, ratio in enumerate(explained_variance_ratio):
    print(f"PC{i+1}: {ratio:.4f}")

# Build a dataframe of PCA coordinates
scores_df = pd.DataFrame(
    scores,
    index=df_t.index,
    columns=[f"PC{i+1}" for i in range(scores.shape[1])]
)

# Define colors for each prefix
prefix_colors = {
    'WT_exp': 'gray',
    'delta5_exp': 'blue',
    'delta9_exp': 'orange',
    'WT_stat': 'black',
    'delta5_stat': 'purple',
    'delta9_stat': 'red'
}

# Plot PCA coordinates of samples
fig, ax = plt.subplots(figsize=(8, 6))

for sample_name in scores_df.index:
    parts = sample_name.split('_')
    prefix = parts[0] + '_' + parts[1]
    color = prefix_colors.get(prefix, 'black')

    ax.scatter(
        scores_df.loc[sample_name, "PC1"],
        scores_df.loc[sample_name, "PC2"],
        color=color
    )


# Add legend
legend_handles = [
    plt.Line2D([0], [0], marker='o', color='w',
               label=prefix, markerfacecolor=color, markersize=10)
    for prefix, color in prefix_colors.items()
]
ax.legend(handles=legend_handles, loc='best', title='Prefix')

ax.set_xlabel(f'PC1 ({explained_variance_ratio[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({explained_variance_ratio[1]*100:.1f}%)')
ax.set_title('PCA of samples')

fig.tight_layout()
fig.savefig('/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/NextSeq/figures/PCA_NAPs_logTPM_plus1_2026.pdf')
plt.show()




### Volcano plot

from bioinfokit import analys, visuz
import pandas as pd
import numpy as np
import os as os

# Change the current working directory
os.chdir('/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/NextSeq/figures')

# load dataset as pandas dataframe
df = pd.read_csv("/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/NextSeq/figures_input/delta9expvsWTexp_IDs.csv")
# Replace empty cells with null values
df.fillna(0, inplace=True)

# list the outliers
outliers = abs(df['log2FoldChange']) >= 6
list_labels = []
list_locus_tags = []
for index, row in df[outliers].iterrows():
    list_labels.append(str(row['Gene Name']))
    list_locus_tags.append(str(row['Locus Tag']))
merged_dict = {}
for i, j in zip(list_locus_tags, list_labels):
    tmp_dict = {i : j}
    if j.startswith('hypothetical'):
        pass
    else:
        merged_dict.update(tmp_dict)

visuz.GeneExpression.volcano(df=df, lfc="log2FoldChange", pv="padj", 
	geneid="Locus Tag", 
    #genenames=(merged_dict),
    show=False,
    gstyle=2, sign_line=True, figname='delta9expvsWTexp_volcano', figtype='pdf', valpha=0.1, 
    #xlm=(-5,5.1,1), ylm=(0,51,10), 
    color=("red", "grey", "blue"),
    lfc_thr=(1, 1), pv_thr=(0.01, 0.01))




### Spike-in analysis

# count all the reads mapped to Bsub and Ecoli separately
cd /Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/NextSeq/sorted_bam_Ecoli
for i in {1..28}; do samtools view -c -F 260 PV$i\_aligned_Ecoli_srt.bam; done
for i in {1..28}; do samtools view -c -F 260 PV$i\_aligned_Bsub_srt.bam; done

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import itertools
import numpy as np

# Load your CSV file
df = pd.read_csv('/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/NextSeq/coverage_sum_up.csv')

# Define a color mapping for the conditions
color_mapping = {
    'WTexp': 'lightgray',
    'WTstat': 'darkgray',
    'delta5exp': 'lightblue',
    'delta5stat': 'darkblue',
    'delta9exp': 'orange',
    'delta9stat': 'darkorange'
}

# Create a new column with the ratio
df['ratio'] = df['number_of_reads_E_coli'] / df['number_of_reads_B_subtilis']

# Extract the condition prefix from the sample names
df['condition'] = df['sample_name'].str.extract(r'(WTexp|WTstat|delta5exp|delta5stat|delta9exp|delta9stat)')

# Specify the desired order for the conditions
order = ['WTexp', 'delta5exp', 'delta9exp', 'WTstat', 'delta5stat', 'delta9stat']
df['condition'] = pd.Categorical(df['condition'], categories=order, ordered=True)

# Group by condition and calculate the mean and standard deviation of the ratios
grouped_data = df.groupby('condition')['ratio'].agg(['mean', 'std', 'count']).loc[order].reset_index()

# Set up the figure and plot the bars
fig, ax = plt.subplots(figsize=(12, 8))
bars = ax.bar(grouped_data['condition'], grouped_data['mean'], 
              yerr=grouped_data['std'], capsize=5, 
              color=[color_mapping[cond] for cond in grouped_data['condition']])
ax.set_ylabel('Mean Ratio (E_coli Reads / B_subtilis Reads)')
#ax.set_title('Mean Ratios by Condition with Standard Deviation')
plt.xticks(rotation=45)

# New function: add a horizontal "stat bar" annotation with the p-value text
def add_bar_annotation(ax, x1, x2, y, text):
    ax.hlines(y, x1, x2, color='k', linewidth=2)
    # Position the text a bit above the line (adjust the multiplier as needed)
    ax.text((x1+x2)/2, y + (0.02 * max_y), text, ha='center', va='bottom')

# Calculate the maximum top of the error bars to help set the base offset
# This is computed as the maximum of (mean + std) for all conditions
max_y = (grouped_data['mean'] + grouped_data['std']).max()
base_offset = max_y * 0.05 if max_y > 0 else 0.1  
used_bars = []  # store (x1, x2, top_y) for previously drawn stat bars

# Define groups to compare: experimental and stationary separately
exp_conditions = ['WTexp', 'delta5exp', 'delta9exp']
stat_conditions = ['WTstat', 'delta5stat', 'delta9stat']

# Function to annotate comparisons for a list of conditions
def annotate_group(conditions):
    global used_bars
    for cond1, cond2 in itertools.combinations(conditions, 2):
        # x positions based on the categorical order
        x1 = order.index(cond1)
        x2 = order.index(cond2)
        
        # Get the ratio data for the two conditions
        data1 = df[df['condition'] == cond1]['ratio']
        data2 = df[df['condition'] == cond2]['ratio']
        
        # Skip if data is missing in either condition
        if len(data1) == 0 or len(data2) == 0:
            continue
        
        # Perform Mann–Whitney U test (non-parametric)
        stat_val, p_val = mannwhitneyu(data1, data2, alternative='two-sided')
        
        # Format the p-value text
        p_text = 'p < 0.001' if p_val < 0.001 else f'p = {p_val:.3f}'
        
        # Calculate the top of the error bars for each condition
        y1 = (grouped_data.loc[grouped_data['condition'] == cond1, 'mean'].values[0] +
              grouped_data.loc[grouped_data['condition'] == cond1, 'std'].values[0])
        y2 = (grouped_data.loc[grouped_data['condition'] == cond2, 'mean'].values[0] +
              grouped_data.loc[grouped_data['condition'] == cond2, 'std'].values[0])
        base = max(y1, y2)
        
        # Determine the y position for the stat bar, adjusting to avoid overlaps
        current_y = base + base_offset
        for (a, b, top_y) in used_bars:
            if not (x2 < a or x1 > b) and current_y <= top_y:
                current_y = top_y + base_offset
        used_bars.append((x1, x2, current_y))
        
        # Draw the horizontal stat bar with p-value text
        add_bar_annotation(ax, x1, x2, current_y, p_text)

# Annotate pairwise comparisons within each group only
annotate_group(exp_conditions)
annotate_group(stat_conditions)

plt.tight_layout()
plt.show()



### Spurious transcription analysis

# Stranded analysis -> only consider reads following the annotations orientations

#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=5:mem=5gb
#PBS -e HTseq_parallel_error.txt

for i in {1..28}
do
   qsub -v index="${i}" /rds/general/user/pvillain/home/RNAseq/NAPs_project/HPC_HTseq_NAPs_RS_stranded.sh
done


#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=5:mem=10gb
#PBS -e HPC_HTseq_NAPs_RS_error.txt

module load anaconda3/personal

source activate /rds/general/user/pvillain/home/anaconda3/envs/HTseq-env

cd /rds/general/user/pvillain/home/RNAseq/NAPs_project

htseq-count -n 1 -t gene -r pos -s no -f bam \
	--idattr=locus_tag \
	--additional-attr=locus_tag \
	--additional-attr=gene \
	--additional-attr=Name \
	--nonunique all \
	--stranded=yes \
	"/rds/general/user/pvillain/home/RNAseq/NAPs_project/Ecoli_aligned/Ecoli_srt_bam/PV${index}_aligned_Ecoli_srt.bam" \
	"/rds/general/user/pvillain/home/RNAseq/NAPs_project/genomes/GCF_015534855.1_genomic.gff" \
	> "/rds/general/user/pvillain/home/RNAseq/NAPs_project/Ecoli_aligned/Ecoli_HTseq_count_stranded/PV${index}_Ecoli_HTseq_stranded.tsv"

htseq-count -n 1 -t gene -r pos -s no -f bam \
	--idattr=locus_tag \
	--additional-attr=locus_tag \
	--additional-attr=gene \
	--additional-attr=Name \
	--nonunique all \
	--stranded=yes \
	"/rds/general/user/pvillain/home/RNAseq/NAPs_project/Bsub_aligned/Bsub_srt_bam/PV${index}_aligned_Bsub_srt.bam" \
	"/rds/general/user/pvillain/home/RNAseq/NAPs_project/genomes/GCF_000009045.1_ASM904v1_genomic.gff" \
	> "/rds/general/user/pvillain/home/RNAseq/NAPs_project/Bsub_aligned/Bsub_HTseq_count_stranded/PV${index}_Bsub_HTseq_stranded.tsv"


# same on reverse annotations

#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=5:mem=5gb
#PBS -e HTseq_parallel_error_reverse.txt

for i in {1..28}
do
   qsub -v index="${i}" /rds/general/user/pvillain/home/RNAseq/NAPs_project/HPC_HTseq_NAPs_RS_stranded_reverse.sh
done


#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=5:mem=10gb
#PBS -e HPC_HTseq_NAPs_RS_error_reverse.txt

module load anaconda3/personal

source activate /rds/general/user/pvillain/home/anaconda3/envs/HTseq-env

cd /rds/general/user/pvillain/home/RNAseq/NAPs_project

htseq-count -n 1 -t gene -r pos -s no -f bam \
	--idattr=locus_tag \
	--additional-attr=locus_tag \
	--additional-attr=gene \
	--additional-attr=Name \
	--nonunique all \
	--stranded=reverse \
	"/rds/general/user/pvillain/home/RNAseq/NAPs_project/Ecoli_aligned/Ecoli_srt_bam/PV${index}_aligned_Ecoli_srt.bam" \
	"/rds/general/user/pvillain/home/RNAseq/NAPs_project/genomes/GCF_015534855.1_genomic.gff" \
	> "/rds/general/user/pvillain/home/RNAseq/NAPs_project/Ecoli_aligned/Ecoli_HTseq_count_stranded_reverse/PV${index}_Ecoli_HTseq_stranded_reverse.tsv"

htseq-count -n 1 -t gene -r pos -s no -f bam \
	--idattr=locus_tag \
	--additional-attr=locus_tag \
	--additional-attr=gene \
	--additional-attr=Name \
	--nonunique all \
	--stranded=reverse \
	"/rds/general/user/pvillain/home/RNAseq/NAPs_project/Bsub_aligned/Bsub_srt_bam/PV${index}_aligned_Bsub_srt.bam" \
	"/rds/general/user/pvillain/home/RNAseq/NAPs_project/genomes/GCF_000009045.1_ASM904v1_genomic.gff" \
	> "/rds/general/user/pvillain/home/RNAseq/NAPs_project/Bsub_aligned/Bsub_HTseq_count_stranded_reverse/PV${index}_Bsub_HTseq_stranded_reverse.tsv"


# Plot the results (anti-sense)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# --- Read and prepare data ---
raw_data = pd.read_csv('/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/Spurious_transcription/HTseq_counts_str_rev_merged_orientation.csv')
raw_data = raw_data.sort_values(by='start')
filtered_data = raw_data[raw_data['convergent'] == 0].copy()

# --- Define conditions and corresponding column names ---
# For each condition, we store the (stranded, reverse) average column names.
conditions_info = {
    'WT_exp': ('WT_exp_stranded_avg', 'WT_exp_reverse_avg'),
    '∆NAP5_exp': ('delta5_exp_stranded_avg', 'delta5_exp_reverse_avg'),
    '∆NAP9_exp': ('delta9_exp_stranded_avg', 'delta9_exp_reverse_avg'),
    'WT_stat': ('WT_stat_stranded_avg', 'WT_stat_reverse_avg'),
    '∆NAP5_stat': ('delta5_stat_stranded_avg', 'delta5_stat_reverse_avg'),
    '∆NAP9_stat': ('delta9_stat_stranded_avg', 'delta9_stat_reverse_avg')
}

# --- Compute aggregated (sum-based) fractions and error bars ---
# We'll compute the central value (p_sum) as:
#    p_sum = (sum(stranded))/(sum(stranded)+sum(reverse))
# and then calculate error bars using the variability from gene-level fractions.
results = {}

for cond, (stranded_col, reverse_col) in conditions_info.items():
    # Aggregated sums and proportion (central value)
    sum_stranded = filtered_data[stranded_col].sum()
    sum_reverse  = filtered_data[reverse_col].sum()
    p_sum = sum_stranded / (sum_stranded + sum_reverse)  # proportion (0–1)
    p_sum_pct = p_sum * 100
    
    # Compute gene-level fractions for each gene
    gene_total = filtered_data[stranded_col] + filtered_data[reverse_col]
    # Avoid division by zero (assumes gene_total > 0 for valid genes)
    p_genes = filtered_data[stranded_col] / gene_total

    # --- Compute error on the arcsine square root scale ---
    # Transform gene-level proportions: t = arcsin(sqrt(p))
    t_genes = np.arcsin(np.sqrt(p_genes))
    sem_t = np.std(t_genes, ddof=1) / np.sqrt(len(t_genes))
    
    # Get the central value on the transformed scale using the aggregated proportion
    t_central = np.arcsin(np.sqrt(p_sum))
    # Determine lower and upper bounds on the transformed scale
    t_lower = t_central - sem_t
    t_upper = t_central + sem_t
    # Back-transform: p = (sin(t))^2
    p_lower = (np.sin(t_lower))**2
    p_upper = (np.sin(t_upper))**2
    # Convert back to percentages
    p_lower_pct = p_lower * 100
    p_upper_pct = p_upper * 100

    # Compute asymmetric error bar lengths (in percentage points)
    err_lower = p_sum_pct - p_lower_pct
    err_upper = p_upper_pct - p_sum_pct

    results[cond] = {
        'p_sum_pct': p_sum_pct,
        'err_lower': err_lower,
        'err_upper': err_upper
    }

# --- Prepare DataFrame for plotting ---
plot_data = pd.DataFrame({
    'Condition': list(results.keys()),
    'AS%': [results[k]['p_sum_pct'] for k in results],
    'err_lower': [results[k]['err_lower'] for k in results],
    'err_upper': [results[k]['err_upper'] for k in results]
})

# Reorder conditions to group by genotype
plot_data['Group'] = ['exp'] * 3 + ['stat'] * 3
order = ['WT_exp', '∆NAP5_exp', '∆NAP9_exp',
         'WT_stat', '∆NAP5_stat', '∆NAP9_stat']

# Unified palette (same color for each genotype, regardless of phase)
shared_palette = {
    'WT_exp': '#CCCCCC',
    '∆NAP5_exp': '#0F80FF',
    '∆NAP9_exp': '#FD8008',
    'WT_stat': '#CCCCCC',
    '∆NAP5_stat': '#0F80FF',
    '∆NAP9_stat': '#FD8008'
}

plt.figure(figsize=(10, 10))
sns.set(style="whitegrid")
ax = sns.barplot(x='Condition', y='AS%', data=plot_data, order=order,
                 palette=[shared_palette[x] for x in order])

# Add asymmetric error bars
for i, row in plot_data.iterrows():
    ax.errorbar(i, row['AS%'], yerr=[[row['err_lower']], [row['err_upper']]], 
                color='black', capsize=5)

# Phase annotation
bar_positions = np.arange(len(order))
# Add horizontal bar and text for 'Exponential'
ax.hlines(y=2, xmin=-0.5, xmax=2.5, color='black', linewidth=1)
ax.text(1, 2.1, 'Exponential', ha='center', va='center', fontsize=16)
# Add horizontal bar and text for 'Stationary'
ax.hlines(y=3.5, xmin=2.5, xmax=5.5, color='black', linewidth=1)
ax.text(4, 3.6, 'Stationary', ha='center', va='center', fontsize=16)

# Aesthetic tweaks
plt.xticks(ticks=bar_positions, labels=order, rotation=45, fontsize=14)
ax.set_ylabel('Proportion of anti-sense transcripts (%)', fontsize=18)
ax.set_xlabel('Condition', fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.savefig("/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/Spurious_transcription/Figures/Anti-sense_transcription_barchart_WT_NAP5_NAP9_expo_stat.pdf", dpi=300, transparent=True)
plt.show()



# Plot the results (intergenic)

import pandas as pd
import matplotlib.pyplot as plt

# ─── Load & prepare ─────────────────────────────────────────────────────────
merged_path = (
    '/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/'
    'Spurious_transcription/intergene_spurious/'
    'Ecoli_HTseq_count_intergene/merged_intergene_scaled_HTseq.csv'
)
df = pd.read_csv(merged_path, index_col=0).apply(pd.to_numeric, errors='coerce')
reads_per_rep = df.sum(axis=0)

# ─── Your mapped replicates per condition ─────────────────────────────────
condition_reps = {
    'WT_exp':      ['WT_exp_2','WT_exp_4','WT_exp_5'],
    'delta5_exp':  ['delta5_exp_2','delta5_exp_3','delta5_exp_4'],
    'delta9_exp':  ['delta9_exp_1','delta9_exp_2','delta9_exp_4'],
    'WT_stat':     ['WT_stat_2','WT_stat_3','WT_stat_5'],
    'delta5_stat': ['delta5_stat_1','delta5_stat_3','delta5_stat_4'],
    'delta9_stat': ['delta9_stat_1','delta9_stat_3','delta9_stat_4']
}

# ─── Compute means & SDs ────────────────────────────────────────────────────
conds, means, stds = [], [], []
for cond, reps in condition_reps.items():
    vals = reads_per_rep[reps]
    conds.append(cond)
    means.append(vals.mean())
    stds.append(vals.std(ddof=1))

# ─── Re‐order & color ───────────────────────────────────────────────────────
desired_order = ['WT_exp','delta5_exp','delta9_exp','WT_stat','delta5_stat','delta9_stat']
idx = [conds.index(c) for c in desired_order]
means = [means[i] for i in idx]
stds  = [stds[i]  for i in idx]
conds = [conds[i] for i in idx]

color_mapping = {
    'WT_exp':     'lightgray',
    'WT_stat':    'darkgray',
    'delta5_exp': 'lightblue',
    'delta5_stat':'blue',
    'delta9_exp': '#FFDAB9',
    'delta9_stat':'orange'
}
colors = [color_mapping[c] for c in conds]

# ─── Plot bars with error bars ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10,7))
x = range(len(conds))
ax.bar(x, means, yerr=stds, capsize=5, color=colors)
ax.set_xticks(x)
ax.set_xticklabels(conds, rotation=90)
ax.set_ylabel('Total intergenic reads (spike-in scaled)')
ax.set_title('Mean ± SD intergenic reads by condition')

# ─── Annotate significant pairs in the STAT group ───────────────────────────
# from your Tukey HSD you found:
#   WT_stat vs delta9_stat  (p≈0.0013)  → '**'
#   delta5_stat vs delta9_stat (p≈0.0009)→ '**'
sig_pairs = [
    ('WT_stat',    'delta9_stat',  '**'),
    ('delta5_stat','delta9_stat',  '**'),
]

# helper to get bar top
bar_tops = [means[i] + stds[i] for i in x]
y_max = max(bar_tops)
y_offset = y_max * 0.05  # 5% of max height

for i, (g1, g2, stars) in enumerate(sig_pairs):
    idx1 = conds.index(g1)
    idx2 = conds.index(g2)
    # horizontal line
    y = y_max + (i+1)*y_offset
    ax.plot([idx1, idx2], [y, y], color='black', linewidth=1.2)
    # small vertical ticks
    ax.plot([idx1, idx1], [y, y-y_offset*0.2], color='black', linewidth=1.2)
    ax.plot([idx2, idx2], [y, y-y_offset*0.2], color='black', linewidth=1.2)
    # stars
    ax.text((idx1+idx2)/2, y + y_offset*0.1, stars, 
            ha='center', va='bottom', fontsize=14)

plt.tight_layout()
plt.savefig("/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/Spurious_transcription/Figures/Intergene_transcription_barchart_WT_NAP5_NAP9_expo_stat.pdf", dpi=300, transparent=True)
plt.show()




### Transcription homogenization


# Scatter plot WT TPM vs log2FC(delta9_exp vs WT_int_exp)

# Load the data
data = pd.read_csv("/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/Transcription_flattening/NAPs_TPM_avg_NaN_DEseq2_WTvsdelta9_exp.csv")

# Convert the 'log2FoldChange' column to numeric, coercing non-numeric values to NaN
data['log2FoldChange'] = pd.to_numeric(data['log2FoldChange'], errors='coerce')

# Drop rows with NaN values in the relevant columns and create a copy
data_clean = data.dropna(subset=['log2FoldChange', 'PV_WT_exp_avg']).copy()

# Define X and Y for the scatter plot
X = data_clean['PV_WT_exp_avg']
Y = data_clean['log2FoldChange']

# Calculate Spearman correlation and p-value
corr, pval = spearmanr(X, Y)

import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# log‐transform X to match the eventual log scale
X = data_clean['PV_WT_exp_avg'].values
Y = data_clean['log2FoldChange'].values
X_log = np.log10(X + 1e-6)   # add a tiny offset to avoid log(0)

# estimate the point density in log-space
xy = np.vstack([X_log, Y])
kde = gaussian_kde(xy, bw_method='scott')   # you can also try bw_method=0.3, etc.
z = kde(xy)

# sort and plot
idx = z.argsort()
X, Y, z = X[idx], Y[idx], z[idx]

plt.figure(figsize=(8, 6))
sc = plt.scatter(
    X, Y,
    c=z, s=10, cmap='magma', alpha=0.7
)
plt.xscale('log')
plt.xlabel("Expression level (TPM)")
plt.ylabel("log2FC of gene expression")
plt.title("")
plt.colorbar(sc, label='Local density')

# annotate correlation
corr, pval = spearmanr(X, Y)
plt.text(
    0.95, 0.95,
    f"Spearman ρ = {corr:.2f}\n p = {pval:.2g}",
    transform=plt.gca().transAxes,
    ha='right', va='top'
)
plt.tight_layout()
plt.savefig("/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/Transcription_flattening/NAP9_scatter_plot_log2FC_vs_WT_TPM.pdf", dpi=300, transparent=True)
plt.show()




### Comparison to the Tjaden transcriptome compendium

import pandas as pd
from scipy.stats import spearmanr

# ==== USER SETTINGS ====
input_file = "/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/Comparison_Tjaden_RNAseq_compendium/Spearman_Paul_2/dataset_dropped_65_TPM.csv"
output_matrix = "/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/Comparison_Tjaden_RNAseq_compendium/Spearman_Paul_2/spearman_correlations_all_vs_all.csv"
output_pvalues = "/Users/pvillain/Documents/Projects/NAPs_project/RNA-seq/Comparison_Tjaden_RNAseq_compendium/Spearman_Paul_2/spearman_pvalues_all_vs_all.csv"
index_col = 0                                # set to None if there's no gene ID column
numeric_only = True                          # drop any non-numeric columns if True

# ==== READ DATA ====
df = pd.read_csv(input_file)
if index_col is not None:
    df = df.set_index(df.columns[index_col])

if numeric_only:
    df = df.select_dtypes(include=["number"])

if df.shape[1] < 2:
    raise SystemExit("Need at least 2 numeric columns to compute correlations.")

# ==== SPEARMAN CORRELATION MATRIX + P-VALUES ====
corr, pvals = spearmanr(df, axis=0, nan_policy="omit")

corr = pd.DataFrame(corr, index=df.columns, columns=df.columns)
pvals = pd.DataFrame(pvals, index=df.columns, columns=df.columns)

# ==== SAVE ====
corr.to_csv(output_matrix)
pvals.to_csv(output_pvalues)

print(f"Saved {corr.shape[0]}x{corr.shape[1]} Spearman correlation matrix to {output_matrix}")
print(f"Saved {pvals.shape[0]}x{pvals.shape[1]} Spearman p-value matrix to {output_pvalues}")


# Plot the correlation to the ∆NAP9 vs the correlation to the WT

from scipy.stats import gaussian_kde
import numpy as np

# Normalize column names once
data_spearman.columns = data_spearman.columns.astype(str).str.strip()
condition_names = data_spearman.columns.tolist()

# Customized scatter plot with specific colors for prefixes and a legend
prefix_colors = {
    'WT_exp': 'blue',
    'WT_stat': 'green',
    'delta5_exp': 'pink',
    'delta9_exp': 'purple',
    'delta5_stat': 'orange',
    'delta9_stat': 'red',
}

# Extract condition names (column headers from filtered_data)
condition_names = data_spearman.columns

# Assign colors based on prefixes
dot_colors = []
for name in condition_names:
    assigned = False
    for prefix, color in prefix_colors.items():
        if name.startswith(prefix):
            dot_colors.append(color)
            assigned = True
            break
    if not assigned:
        dot_colors.append('black')

# Style (match your reference)
plt.style.use("default")
plt.rcParams.update({
    "axes.grid": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

plt.figure(figsize=(8, 6))

# Split indices
black_idx = [i for i, c in enumerate(dot_colors) if c == 'black']
color_idx = [i for i, c in enumerate(dot_colors) if c != 'black']

# === KDE density for black points ===
x_black = data_spearman['WT_exp_1'].iloc[black_idx].values
y_black = data_spearman['delta9_exp_1'].iloc[black_idx].values

xy = np.vstack([x_black, y_black])
z = gaussian_kde(xy)(xy)

# sort so low density plotted first
idx = z.argsort()
x_black, y_black, z = x_black[idx], y_black[idx], z[idx]

# Plot density-colored background
sc = plt.scatter(
    x_black,
    y_black,
    c=z,
    s=10,
    cmap='magma',
    alpha=0.7
)

plt.colorbar(sc, label='Local density')

# === Plot highlighted colored points LAST ===
plt.scatter(
    data_spearman['WT_exp_1'].iloc[color_idx],
    data_spearman['delta9_exp_1'].iloc[color_idx],
    c=[dot_colors[i] for i in color_idx],
    alpha=0.9,
    s=20,
    zorder=3
)

# Legend
for prefix, color in prefix_colors.items():
    plt.scatter([], [], color=color, label=prefix, alpha=0.9, s=50)

plt.legend(
    title="Color Key (Prefixes)",
    title_fontsize=12,
    fontsize=10,
    loc='upper left',
    frameon=True
)

plt.title(
    'Spearman rho to delta9_exp_1 vs Spearman rho to WT_exp_1',
    fontsize=16
)
plt.xlabel('Spearman rho to WT_exp_1', fontsize=14)
plt.ylabel('Spearman rho to delta9_exp_1', fontsize=14)

# ensure no gridlines
ax = plt.gca()
ax.grid(False, which='both', axis='both')

plt.tight_layout()
plt.savefig("/Users/pvillain/Documents/Manuscripts/NAPs_paper/Figures_transcription/Figures_comparison_Tjaden_dataset/Figure_scatter_rho_delta9_exp_1_vs_WT_exp_1.pdf", dpi=300, transparent=True)
plt.show()


