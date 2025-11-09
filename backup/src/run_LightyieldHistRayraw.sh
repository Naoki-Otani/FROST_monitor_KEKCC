#!/bin/bash

# 実行対象のファイル名（引数）を配列で定義
file_list=(
  # "run00007"
  "run00010"
  # "run00012"
  # "run00014"
  # "run00016"
  # "run00018"
  "run00020"
  "run00022"
  "run00030"
  "run00032"
  "run00034"
)

# 実行ループ
for fname in "${file_list[@]}"; do
  echo "Running LightyieldHistRayraw for ${fname} ..."
  root -l -b -q "LightyieldHistRayraw.C(\"${fname}\")"
done

echo "All LightyieldHistRayraw runs completed."
