#!/bin/bash

# 実行対象のファイル名（引数）を配列で定義
file_list=(
    #"run00004"
    # "run00006"
    # "run00009"
    # "run00011"
    # "run00013"
    # "run00015"
    # "run00017"
    # "run00019"
    # "run00021"
    "run00029"
    "run00031"
    "run00033"

)

# 実行ループ
for fname in "${file_list[@]}"; do
  echo "Running calibrayraw for ${fname} ..."
  #root -l -b -q "calibrayraw.C(\"${fname}\",\"chmap_run40\")"
  root -l -b -q "calibrayraw.C(\"${fname}\",\"chmap_20251009\")"
done

echo "All calibrayraw runs completed."
