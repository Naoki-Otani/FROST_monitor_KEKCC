#!/bin/bash

# 実行対象のファイル名（引数）を配列で定義
file_list=(
  "run00044"
  "run00046"
  "run00048"
  "run00050"
  "run00051"
  "run00056"
  "run00058"
  "run00061"
  "run00062"
  "run00071"
)

# 実行ループ
for fname in "${file_list[@]}"; do
  echo "Running adcint_rayraw for ${fname} ..."
  root -l -b -q "adcint_rayraw.C(\"${fname}\",\"chmap_run40\")"
done

echo "All adcint_rayraw runs completed."
