#!/bin/bash

# 実行対象のファイル名（引数）を配列で定義
file_list=(
  # "run44"
  # "20250923-1"
  # "run46"
  # "20250923-2"
  # "run48"
  # "20250924-1"
  # "run50"
  # "20250924-2"
  # "run51"
  # "20250924-2"
  # "run56"
  # "20250924-2"
  # "run58"
  # "20250925-1"
  # "run61"
  # "20250925-1"
  # "run62"
  # "20250925-1"
  "run63"
  "20250925-1"
  "run64"
  "20250925-1"
  "run67"
  "20250926-1"
  "run68"
  "20250926-1"
  "run69"
  "20250926-1"
  # "run71"
  # "20250928-1"
  "run72"
  "20250928-1"
)

# 2つずつペアで実行
for ((i=0; i<${#file_list[@]}; i+=2)); do
  fname1=${file_list[$i]}
  fname2=${file_list[$i+1]}

  echo "Running convertlightyield_caen for: $fname1 and $fname2"
  root -l -b -q "convertlightyield_caen.C(\"$fname1\", \"$fname2\")"
done

echo "✅ All convertlightyield_caen.C runs completed."
