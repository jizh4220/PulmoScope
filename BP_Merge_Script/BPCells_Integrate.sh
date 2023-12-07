#!/bin/bash
rds_list_path=$1
merged_h5ad_path=$2

echo "Running Bulk BPCells Integration"
source activate BPCells;

echo "Running rds2h5ad.R..."
Rscript rds2h5ad.R $rds_list_path;

echo "Running scanpyMergeCounts.py..."
python3 scanpyMergeCounts.py $rds_list_path $merged_h5ad_path F ;
# Once Scanpy h5ad file is successfully prepared, call BPCell_Preparation
echo "Running BPCell_Preparation.R..."
Rscript BPCell_Preparation.R $merged_h5ad_path;

echo "Integrated BPCells object has been prepared successfully!"
