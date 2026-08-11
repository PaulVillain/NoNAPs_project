### Reads processing

#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=8:mem=20gb
#PBS -e trimm_error.txt

eval "$(~/anaconda3/bin/conda shell.bash hook)"

source activate /rds/general/user/pvillain/home/anaconda3/envs/TrimGalore-env

cd /rds/general/user/pvillain/home/NAPs_genomes/Marker_frequency_analysis/ONT_and_Illumina

trim_galore --cores 8 --paired MQ4KBQ_1_WT_illumina_R1.fastq MQ4KBQ_1_WT_illumina_R2.fastq
trim_galore --cores 8 --paired MQ4KBQ_2_delta5_illumina_R1.fastq MQ4KBQ_2_delta5_illumina_R2.fastq
trim_galore --cores 8 --paired MQ4KBQ_3_delta9_illumina_R1.fastq MQ4KBQ_3_delta9_illumina_R2.fastq


qsub HPC_TrimGalore.sh


eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate /rds/general/user/pvillain/home/anaconda3/envs/bowtie2
bowtie2-build MQ4KBQ_1_WT_1_reference.fna WT
bowtie2-build MQ4KBQ_2_delta5_1_reference.fna delta5
bowtie2-build MQ4KBQ_3_delta9_1_reference.fna delta9
    

#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=30:mem=50gb
#PBS -e bowtie2_MFA_error.txt

eval "$(~/anaconda3/bin/conda shell.bash hook)"

source activate /rds/general/user/pvillain/home/anaconda3/envs/bowtie2

cd /rds/general/user/pvillain/home/NAPs_genomes/Marker_frequency_analysis/ONT_and_Illumina

bowtie2 -x /rds/general/user/pvillain/home/NAPs_genomes/Marker_frequency_analysis/ONT_and_Illumina/WT -1 MQ4KBQ_1_WT_illumina_R1_val_1.fq -2 MQ4KBQ_1_WT_illumina_R2_val_2.fq -S /rds/general/user/pvillain/home/NAPs_genomes/Marker_frequency_analysis/ONT_and_Illumina/sam_files/WT_illumina.sam
bowtie2 -x /rds/general/user/pvillain/home/NAPs_genomes/Marker_frequency_analysis/ONT_and_Illumina/delta5 -1 MQ4KBQ_2_delta5_illumina_R1_val_1.fq -2 MQ4KBQ_2_delta5_illumina_R2_val_2.fq -S /rds/general/user/pvillain/home/NAPs_genomes/Marker_frequency_analysis/ONT_and_Illumina/sam_files/delta5_illumina.sam
bowtie2 -x /rds/general/user/pvillain/home/NAPs_genomes/Marker_frequency_analysis/ONT_and_Illumina/delta9 -1 MQ4KBQ_3_delta9_illumina_R1_val_1.fq -2 MQ4KBQ_3_delta9_illumina_R2_val_2.fq -S /rds/general/user/pvillain/home/NAPs_genomes/Marker_frequency_analysis/ONT_and_Illumina/sam_files/delta9_illumina.sam

qsub HPC_Bowtie2_MFA.sh



#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=10:mem=50gb
#PBS -e samtools_reads_error.txt


eval "$(~/anaconda3/bin/conda shell.bash hook)"

source activate /rds/general/user/pvillain/home/anaconda3/envs/samtools-env

cd /rds/general/user/pvillain/home/NAPs_genomes/Marker_frequency_analysis/ONT_and_Illumina/sam_files

samtools sort -@ 10 WT_illumina.sam -o WT_illumina.bam
samtools sort -@ 10 delta5_illumina.sam -o delta5_illumina.bam
samtools sort -@ 10 delta9_illumina.sam -o delta9_illumina.bam

samtools sort -@ 10 WT_illumina.bam -o WT_illumina_srt.bam
samtools sort -@ 10 delta5_illumina.bam -o delta5_illumina_srt.bam
samtools sort -@ 10 delta9_illumina.bam -o delta9_illumina_srt.bam

samtools index -@ 10 WT_illumina_srt.bam
samtools index -@ 10 delta5_illumina_srt.bam
samtools index -@ 10 delta9_illumina_srt.bam


qsub HPC_samtools_MFA_illumina.sh
    


conda activate /Users/pvillain/miniconda3/envs/Deeptools_env
    
bamCoverage --binSize 1 --normalizeUsing BPM -b WT_illumina_srt.bam -o WT_illumina_srt.bw
bamCoverage --binSize 1 --normalizeUsing BPM -b delta5_illumina_srt.bam -o delta5_illumina_srt.bw
bamCoverage --binSize 1 --normalizeUsing BPM -b delta9_illumina_srt.bam -o delta9_illumina_srt.bw


bamCoverage --binSize 1 --normalizeUsing BPM -b aln_MQ4KBQ_1_WT_nanopore.srt.bam -o WT_ONT_srt.bw
bamCoverage --binSize 1 --normalizeUsing BPM -b aln_MQ4KBQ_1_delta5_nanopore.srt.bam -o delta5_ONT_srt.bw
bamCoverage --binSize 1 --normalizeUsing BPM -b aln_MQ4KBQ_1_delta9_nanopore.srt.bam -o delta9_ONT_srt.bw





### Normalize the coverage

import pyBigWig
import numpy as np
import os

os.chdir("/Users/pvillain/Documents/Projects/NAPs_project/Marker_frequency_analysis/Plasmidsaurus_hybrid_sequencing/ONT_and_illumina")

bw = pyBigWig.open("WT_illumina_srt.bw")

chrom = list(bw.chroms().keys())[0]   # contig_1
L = bw.chroms()[chrom]

bin_size = 1000
bin_means = []

for start in range(0, L, bin_size):
    end = min(start + bin_size, L)
    v = bw.values(chrom, start, end)
    bin_means.append(np.nanmean(v))

genome_mean_WT = np.nanmean(bin_means)
mfa_WT = np.array(bin_means) / genome_mean_WT

print("chrom:", chrom, "length:", L, "genome_mean_WT:", genome_mean_WT)
print("first 5 MFA:", mfa_WT[:5])

bw = pyBigWig.open("delta5_illumina_srt.bw")

chrom = list(bw.chroms().keys())[0]   # contig_1
L = bw.chroms()[chrom]

bin_size = 1000
bin_means = []

for start in range(0, L, bin_size):
    end = min(start + bin_size, L)
    v = bw.values(chrom, start, end)
    bin_means.append(np.nanmean(v))

genome_mean_NAP5 = np.nanmean(bin_means)
mfa_NAP5 = np.array(bin_means) / genome_mean_NAP5

print("chrom:", chrom, "length:", L, "genome_mean_NAP5:", genome_mean_NAP5)
print("first 5 MFA:", mfa_NAP5[:5])

bw = pyBigWig.open("delta9_illumina_srt.bw")

chrom = list(bw.chroms().keys())[0]   # contig_1
L = bw.chroms()[chrom]

bin_size = 1000
bin_means = []

for start in range(0, L, bin_size):
    end = min(start + bin_size, L)
    v = bw.values(chrom, start, end)
    bin_means.append(np.nanmean(v))

genome_mean_NAP9 = np.nanmean(bin_means)
mfa_NAP9 = np.array(bin_means) / genome_mean_NAP9

print("chrom:", chrom, "length:", L, "genome_mean_NAP9:", genome_mean_NAP9)
print("first 5 MFA:", mfa_NAP9[:5])



### Figure with just the Illumina reads

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === Fill these in (0-based bp coordinates on each contig) ===
ori_WT   = 4589927      # <-- replace
ori_NAP5 = 4588716      # <-- replace
ori_NAP9 = 4587244      # <-- replace

bin_size = 1000
L_WT = 4632016
L_5  = 4630807
L_9  = 4629334

def build_df(mfa, genome_length, label):
    starts = np.arange(0, genome_length, bin_size)
    mids = starts + bin_size / 2
    mids[mids > genome_length] = genome_length

    return pd.DataFrame({
        "mid": mids,
        "mfa": np.asarray(mfa, dtype=float),
        "strain": label
    })

def smooth(y, window=100):
    return (
        pd.Series(y)
        .rolling(window=window, center=True, min_periods=1)
        .mean()
        .to_numpy()
    )

def rotate_to_ori(df, genome_length, ori_pos):
    """
    Rotate circular genome so that ori_pos becomes 0, then recenter to [-L/2, +L/2).
    Returns x in Mb.
    """
    x = df["mid"].to_numpy()

    # rotate so ori is at 0 in [0, L)
    x_rot = (x - ori_pos) % genome_length

    # recenter to [-L/2, +L/2)
    x_center = np.where(x_rot > genome_length/2, x_rot - genome_length, x_rot)

    df2 = df.copy()
    df2["x_mb"] = x_center / 1e6
    return df2.sort_values("x_mb")

def rotate_coord_to_ori_mb(coord, genome_length, ori_pos):
    """
    Transform a single genomic coordinate (bp, original reference) into the
    rotated/recentered coordinate system used for plotting (Mb).
    """
    x_rot = (coord - ori_pos) % genome_length
    if x_rot > genome_length / 2:
        x_rot -= genome_length
    return x_rot / 1e6

# Build per-strain dfs (uses your already computed mfa_* arrays and L_* lengths)
df_WT   = build_df(mfa_WT,   L_WT, "WT")
df_NAP5 = build_df(mfa_NAP5, L_5,  "ΔNAP5")
df_NAP9 = build_df(mfa_NAP9, L_9,  "ΔNAP9")

# Rotate to oriC
df_WT_r   = rotate_to_ori(df_WT,   L_WT, ori_WT)
df_NAP5_r = rotate_to_ori(df_NAP5, L_5,  ori_NAP5)
df_NAP9_r = rotate_to_ori(df_NAP9, L_9,  ori_NAP9)

# Plot separately (3 panels)
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(4, 8), sharey=True)

axes[0].scatter(df_WT_r["x_mb"], df_WT_r["mfa"], s=3, alpha=0.3)
axes[0].plot(df_WT_r["x_mb"], smooth(df_WT_r["mfa"]), color="black", linewidth=1)
axes[0].set_title("WT")

axes[1].scatter(df_NAP5_r["x_mb"], df_NAP5_r["mfa"], s=3, alpha=0.3)
axes[1].plot(df_NAP5_r["x_mb"], smooth(df_NAP5_r["mfa"]), color="black", linewidth=1)
axes[1].set_title("ΔNAP5")

axes[2].scatter(df_NAP9_r["x_mb"], df_NAP9_r["mfa"], s=3, alpha=0.3)
axes[2].plot(df_NAP9_r["x_mb"], smooth(df_NAP9_r["mfa"]), color="black", linewidth=1)
axes[2].set_title("ΔNAP9")

# ---- Highlight interval on ΔNAP9 using ORIGINAL (pre-rotation) coordinates ----
start_raw, end_raw = 2441369, 2470146  # original bp coordinates (0-based)
x_start = rotate_coord_to_ori_mb(start_raw, L_9, ori_NAP9)
x_end   = rotate_coord_to_ori_mb(end_raw,   L_9, ori_NAP9)
xmin, xmax = sorted([x_start, x_end])

y_highlight = 1.4  # pick the MFA level where you want the horizontal line
axes[2].hlines(
    y=y_highlight,
    xmin=xmin,
    xmax=xmax,
    colors="red",
    linewidth=2
)


for ax in axes:
    ax.axhline(1, linestyle="--", linewidth=1, alpha=0.5)
    ax.set_ylabel("Normalized read count")

axes[-1].set_xlabel("Genome position relative to oriC (Mb)")

for ax in axes:
    ax.set_ylim(0, 2.0)

plt.tight_layout()
plt.savefig("/Users/pvillain/Documents/Projects/NAPs_project/Marker_frequency_analysis/Figures/MFA_WT_NAP5_NAP9_ONT_only_300dpi.pdf", dpi=300, transparent=True)
plt.show()
