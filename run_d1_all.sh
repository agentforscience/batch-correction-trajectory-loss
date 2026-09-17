#!/bin/bash
cd /workspaces/how_much_biology_do_batch_corr_20260906_074942_ba309bd5
source .venv/bin/activate
export OMP_NUM_THREADS=6 MKL_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 NUMEXPR_NUM_THREADS=6
for i in 0 1 2 3; do
  gpu=$(( i % 2 == 0 ? 1 : 3 ))
  CUDA_VISIBLE_DEVICES=$gpu python src/run_d1.py --seeds 1 2 3 4 5 --slice $i 4 \
     --out results/d1_raw_w$i.csv > logs/d1_w$i.log 2>&1 &
done
wait
echo "ALL_D1_DONE"
