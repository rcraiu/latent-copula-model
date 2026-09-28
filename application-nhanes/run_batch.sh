#!/bin/bash
# Launch several chains in parallel and block until all have finished.
#
# Run from the repository root:
#     M=600000 ./application-nhanes/run_batch.sh Gaussian:1 Gaussian:2 Frank:1
#
# Each R process is single-threaded, so that running one per family and chain
# does not oversubscribe the cores. Ten chains on eight cores is a reasonable
# load; beyond that, split into batches.

set -u

# resolve the repository root from this script's location, then work from there
cd "$(dirname "$0")/.." || exit 1

if [ "$#" -eq 0 ]; then
  echo "usage: [M=600000] $0 <Family>:<chain> [<Family>:<chain> ...]" >&2
  echo "families: Gaussian Clayton Frank Gumbel Independence" >&2
  exit 2
fi

export OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1

M=${M:-600000}
ENGINE=${ENGINE:-cpp}
PROPOSAL=${PROPOSAL:-block}
BETA=${BETA:-norm1}

mkdir -p output/logs output/fits

echo "M=$M engine=$ENGINE proposal=$PROPOSAL beta=$BETA"
for spec in "$@"; do
  fam=${spec%%:*}
  ch=${spec##*:}
  Rscript application-nhanes/fit_nhanes.R "$fam" "$ch" "$M" "$ENGINE" "$PROPOSAL" "$BETA" \
      > "output/logs/${fam}_${ch}.log" 2>&1 &
  echo "launched $fam chain $ch as pid $!"
done

wait
echo "BATCH COMPLETE"
ls -la output/fits/
