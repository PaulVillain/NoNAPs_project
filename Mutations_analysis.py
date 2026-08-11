### Mapping of the read on the reference genome (minimap2) - same process for the ∆NAP strains and the evolved lineages

#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=10:mem=12gb
#PBS -e Minimap2_LTC50_error.txt
#PBS -o Minimap2_LTC50_output.txt

module load anaconda3/personal

source activate /rds/general/user/pvillain/home/anaconda3/envs/nanopore_seq_env

cd /rds/general/user/pvillain/home/NAPs_genomes/LTC50

# Map the read on the ref genome

for i in {1..5}; do minimap2 -c -ax map-ont /rds/general/user/pvillain/home/NAPs_genomes/LTC50/NZ_CP064677.fasta /rds/general/user/pvillain/home/NAPs_genomes/LTC50/LTC50_WT_$i\.fastq -t 10 > /rds/general/user/pvillain/home/NAPs_genomes/LTC50/sam_files/aln_LTC50_WT_$i\.sam; done
for i in {1..5}; do minimap2 -c -ax map-ont /rds/general/user/pvillain/home/NAPs_genomes/LTC50/NZ_CP064677.fasta /rds/general/user/pvillain/home/NAPs_genomes/LTC50/LTC50_NAP9_$i\.fastq -t 10 > /rds/general/user/pvillain/home/NAPs_genomes/LTC50/sam_files/aln_LTC50_NAP9_$i\.sam; done

qsub HPC_minimap2_LTC50.sh




### Index and sort the sam files (SamTools)

conda activate /Users/pvillain/ENTER/envs/samtools-env

# SAM TO BAM conversion
for i in {1..5}; do samtools view -@ 6 -bS aln_LTC50_WT_$i\.sam > aln_LTC50_WT_$i\.bam; done
for i in {1..5}; do samtools view -@ 6 -bS aln_LTC50_NAP9_$i\.sam > aln_LTC50_NAP9_$i\.bam; done

# Sort your bam file
for i in {1..5}; do samtools sort -@ 6 aln_LTC50_WT_$i\.bam -o aln_LTC50_WT_$i\_srt.bam; done
for i in {1..5}; do samtools sort -@ 6 aln_LTC50_NAP9_$i\.bam -o aln_LTC50_NAP9_$i\_srt.bam; done

# Index your bam file
for i in {1..5}; do samtools index aln_LTC50_WT_$i\_srt.bam; done
for i in {1..5}; do samtools index aln_LTC50_NAP9_$i\_srt.bam; done

conda activate /Users/pvillain/miniconda3/envs/Deeptools_env
# get the TPM coverage for IGV visualisation
for i in {1..5}; do bamCoverage --binSize 1 --normalizeUsing BPM -b aln_LTC50_WT_$i\_srt.bam -o LTC50_WT_$i\_srt.bw; done
for i in {1..5}; do bamCoverage --binSize 1 --normalizeUsing BPM -b aln_LTC50_NAP9_$i\_srt.bam -o LTC50_NAP9_$i\_srt.bw; done




### Mutation calling  (Clair3)

conda activate /Users/pvillain/ENTER/envs/samtools-env

samtools faidx /Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/NZ_CP064677.fasta

### Warning: the fasta reference has to be in the same folder than the bam files, otherwise the software crash!!! ###

for i in {1..5}; do docker run -it \
  -v /Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/bam_files:/Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/bam_files \
  -v /Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/Clair3_output_haploid_sensitive:/Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/Clair3_output_haploid_sensitive \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=/Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/bam_files/aln_LTC50_WT_$i\_srt.bam \
  --ref_fn=/Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/bam_files/NZ_CP064677.fasta \
  --threads=4 \
  --platform="ont" \
  --model_path="/opt/models/r941_prom_sup_g5014" \
  --output=/Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/Clair3_output_haploid_sensitive/WT_$i\_Clair3 \
  --no_phasing_for_fa \
  --include_all_ctgs \
  --enable_long_indel \
  --haploid_sensitive; done

for i in {1..5}; do docker run -it \
  -v /Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/bam_files:/Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/bam_files \
  -v /Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/Clair3_output_haploid_sensitive:/Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/Clair3_output_haploid_sensitive \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=/Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/bam_files/aln_LTC50_NAP9_$i\_srt.bam \
  --ref_fn=/Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/bam_files/NZ_CP064677.fasta \
  --threads=4 \
  --platform="ont" \
  --model_path="/opt/models/r941_prom_sup_g5014" \
  --output=/Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/Clair3_output_haploid_sensitive/NAP9_$i\_Clair3 \
  --no_phasing_for_fa \
  --include_all_ctgs \
  --enable_long_indel \
  --haploid_sensitive; done




### Intersect Clair3 output with the gff information

conda activate bedtools-env

# get the mutation information + the gff information for the overlaps
for i in {1..5}; do bedtools intersect -header -wo -a /Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/bam_files/GCF_015534855.1_genomic_NZ_CP064677.gff -b /Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/Clair3_output/merge_output_WT_$i\.vcf > /Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/intersect_gff_Clair3/intersect_gff_wo_WT_$i\.txt; done
for i in {1..5}; do bedtools intersect -header -wo -a /Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/bam_files/GCF_015534855.1_genomic_NZ_CP064677.gff -b /Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/Clair3_output/merge_output_NAP9_$i\.vcf > /Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/intersect_gff_Clair3/intersect_gff_wo_NAP9_$i\.txt; done


# Clean the output
import csv

# 1) Adjust these paths as needed:
in_path  = "/Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/intersect_gff_Clair3/intersect_gff_WT_5.txt"
out_path = "/Users/pvillain/Documents/Projects/NAPs_project/LTC50_experiment/Plasmidsaurus_sequencing/intersect_gff_Clair3/intersect_gff_WT_5_clean.tsv"

# 2) Names to exclude:
EXCLUDE = {'hupA','hupB','hns','stpA','fis','lrp','dps','ihfA','ihfB'}

def parse_attrs(attr_str):
    """Parse the 9th GFF column into a dict of {key: value}."""
    d = {}
    for part in attr_str.split(';'):
        if '=' in part:
            k, v = part.split('=', 1)
            d[k] = v
    return d

with open(in_path, newline='') as in_fh, open(out_path, 'w', newline='') as out_fh:
    reader = csv.reader(in_fh, delimiter='\t')
    writer = csv.writer(out_fh, delimiter='\t')

    for row in reader:
        # pass through comments or blank lines
        if not row or row[0].startswith('#') or len(row) < 9:
            writer.writerow(row)
            continue

        attrs = parse_attrs(row[8])
        gene = attrs.get('gene', '')
        name = attrs.get('Name', '')

        # only write rows where neither gene nor Name is in EXCLUDE
        if gene not in EXCLUDE and name not in EXCLUDE:
            writer.writerow(row)

print(f"Done → filtered file written to “{out_path}”")
