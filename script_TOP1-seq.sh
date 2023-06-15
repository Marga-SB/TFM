#!/bin/bash

#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -t 99:00:00
#SBATCH -o log_wes.out
#SBATCH -e log_wes.err
#SBATCH -p ladon
#SBATCH -J wes_baranello

export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_EXECUTOR=slurm
# Set dir to store singularity images
export NXF_SINGULARITY_CACHEDIR="/mnt/beegfs/home/marsabbon/.nxf_singularity_cache"
#ulimit -u 4126507
#ulimit -c unlimited

# Get slurm requested RAM and set MEM env variable
#MEM=${SLURM_MEM_PER_NODE}


# Set tmp dirs to prevent /tmp from being filled up
export TMPDIR=/scratch/marsabbon/singularity/tmp/$SLURM_JOB_ID
export SINGULARITY_TMPDIR=/scratch/marsabbon/singularity/tmp/$SLURM_JOB_ID

mkdir -p $TMPDIR

# Copy files to /scratch partition on nodes or create symlinks otherwise
SCRATCH=true
#SCRATCH=false
STAGEINMODE="copy"
#STAGEINMODE="symlink"

# Launch nextflow
## Chipseq pipeline
## Genomes are stored in /mnt/beegfs/genomes/igenomes/
nextflow -Dnxf.pool.type=sync run /mnt/beegfs/home/marsabbon/.nxf_singularity_cache/nf-core-chipseq-2.0.0/workflow -w work -process.errorStrategy='ignore' -process.scratch=${SCRATCH} -process.stageInMode=${STAGEINMODE} -process.maxForks=80 -process.time='2h' -process.cache='lenient' -profile singularity --input samples/samplesheet_input/baranello.csv --genome "GRCh38" --fasta /mnt/beegfs/genomes/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/version0.6.0/genome.fa --gtf /mnt/beegfs/genomes/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/genes.gtf --outdir results_baranello_broad_peak --broad --read_length 50 --macs_fdr 0.05  
