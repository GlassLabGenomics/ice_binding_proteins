#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --partition=t1small
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=100MB
#SBATCH --time=03:00:00
#SBATCH --array=1-2

set -o errexit
#set -o nounset
module --quiet purge 

module load GCC/11.3.0
module load OpenMPI/4.1.4
module load BLAST+/2.13.0

### Read in fastas
INFILE=filelist
findfasta=$(awk "NR==$SLURM_ARRAY_TASK_ID" $INFILE)

### Set I/O paths
QUERY="fastas/${findfasta}"
XMLOUT="blastxmls/${findfasta::-6}.xml"

echo "Input fasta: ${QUERY}"
echo "Output XML written to: ${XMLOUT}"

echo "time blastp -query ${QUERY} -db /center1/GLASSLAB/yhsieh/blastdb/refseq_protein -out ${XMLOUT} -evalue 0.001 -outfmt 16 -num_threads 16 -max_target_seqs 5000"

time blastp -query ${QUERY} -db /center1/GLASSLAB/yhsieh/blastdb/refseq_protein -out ${XMLOUT} -evalue 0.001 -outfmt 16 -num_threads 16 -max_target_seqs 5000
