#!/bin/bash

# 対象のファイル一覧
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

# weight の値をここに定義（何個でもOK）
weight_list=(
  1.0
  1.5
  2.0
  2.5
  3.0
  3.5
  4.0
)

# 二重ループで両方回す
for fname in "${file_list[@]}"; do
  for weight in "${weight_list[@]}"; do
    echo "Running xgyg for ${fname} with weight=${weight} ..."
    root -l -b -q "xgyg.C(\"${fname}\", ${weight})"
  done
done

echo "All xgyg runs completed."
