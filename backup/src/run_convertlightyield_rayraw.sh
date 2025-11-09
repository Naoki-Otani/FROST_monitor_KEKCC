#!/bin/bash

# 実行対象のファイル名（引数）を配列で定義
file_list=(
  # "run00007"
  # "run00006"
  "run00010"
  "run00009"
  # "run00012"
  # "run00011"
  # "run00014"
  # "run00013"
  # "run00016"
  # "run00015"
  # "run00018"
  # "run00017"
  "run00020"
  "run00019"
  "run00022"
  "run00021"
  "run00030"
  "run00029"
  "run00032"
  "run00031"
  "run00034"
  "run00033"
)

# 2つずつペアで実行
for ((i=0; i<${#file_list[@]}; i+=2)); do
  fname1=${file_list[$i]}
  fname2=${file_list[$i+1]}

  echo "Running convertlightyield_rayraw for: $fname1 and $fname2"
  root -l -b -q "convertlightyield_rayraw.C(\"$fname1\",\"chmap_20251009\",\"ReferenceGain_fiberdif\", \"$fname2\")"
done

echo "✅ All convertlightyield_rayraw.C runs completed."
