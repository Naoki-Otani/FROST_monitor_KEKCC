#!/bin/bash

# 実行対象のファイル名（引数）を配列で定義
file_list=(
  "20251009-1"
  "20251010-1"
  "20251011-1"
  "20251011-2"
  "20251012-1"
  "20251012-2"
  "20251013-1"
)

# 実行ループ
for fname in "${file_list[@]}"; do
  echo "Running calibcaen for ${fname} ..."
  root -l -b -q "calibcaen.C(\"${fname}\")"
done

echo "All calibcaen runs completed."
