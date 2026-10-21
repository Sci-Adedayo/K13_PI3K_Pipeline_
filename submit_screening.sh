#!/bin/bash
mkdir -p logs

NUM_MUTANTS=$(ls -1 mutants/*_A.pdb 2>/dev/null | wc -l)

if [ "$NUM_MUTANTS" -eq 0 ]; then
    echo "ERROR: No mutant structures found in mutants/*_A.pdb"
    exit 1
fi

TOTAL_TASKS=$(( NUM_MUTANTS * 2 ))
MAX_ARRAY_INDEX=$(( TOTAL_TASKS - 1 ))

echo "Found $NUM_MUTANTS mutants across 2 domains ($TOTAL_TASKS total docking runs)."
echo "Submitting SLURM array: 0-$MAX_ARRAY_INDEX"

sbatch --array=0-${MAX_ARRAY_INDEX} run_dual_domain_array.slurm
