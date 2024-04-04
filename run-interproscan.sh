#!/bin/bash
#SBATCH --mem=64G
#SBATCH --output=interpro_results.out
#SBATCH --mail-user=bdm50@duke.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

cd /work/bdm50/new_interpro/interproscan-5.65-97.0/lib/

conda info --env

./interproscan.sh -i ./Hery_peptide_models.fasta -goterms -cpu 16