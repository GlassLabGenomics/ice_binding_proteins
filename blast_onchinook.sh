#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --partition=t1small
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=100MB
#SBATCH --time=03:00:00

set -o errexit
set -o nounset

module --quiet purge 
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load BLAST+/2.13.0

QUERY=$1
XMLOUT=$2

time blastp -query ${QUERY} -db /center1/GLASSLAB/yhsieh/blastdb/refseq_protein -out ${XMLOUT} -evalue 0.001 -outfmt 16 -num_threads 16 -max_target_seqs 5000
