#!/bin/bash
set -euo pipefail

SCRIPT=apply_selection.py
NPROCS=4   # use 2–4 on lxplus

YEARS=(2016 2017 2018)
MAGS=(up down)
COMPS=(
  b0_dstmunu
  b0_dsttau
  b0_dststtau
  b+_dststtau
  b+_dstdx
  b0_dstdx
  b0_dststmu_415
  b0_dststmu_10413
  b0_dststmu_20413
  b+_dststmu_425
  b+_dststmu_10423
  b+_dststmu_20423
  bs_dststmu_10433
  bs_dststmu_435
  b+_dststmu_heavy
  b0_dststmu_heavy
)

mkdir -p logs/sels

for YEAR in "${YEARS[@]}"; do
  for MAG in "${MAGS[@]}"; do
    for COMP in "${COMPS[@]}"; do

      LOG="logs/sels/${YEAR}_${MAG}_${COMP}.log"
      echo "[START] ${YEAR} ${MAG} ${COMP}"

      python3 "${SCRIPT}" "${YEAR}" "${MAG}" "${COMP}" \
        > "$LOG" 2>&1 &

      # HARD throttle: wait until a slot is free
      while (( $(jobs -r | wc -l) >= NPROCS )); do
        sleep 1
      done

    done
  done
done

wait
echo "✅ All jobs finished"
