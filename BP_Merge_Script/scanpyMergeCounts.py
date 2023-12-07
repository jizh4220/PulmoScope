import pandas as pd
import scanpy as sc
import scipy
import glob
import os
import sys
import re
import anndata
import numpy as np
# import datatable as dt

def merge_duplicate_gene(dupl_gene_counts, rest_counts):
    """
    Merges the expression counts for duplicated genes.

    :param dupl_gene_list: List of gene names that are duplicated.
    :param dupl_gene_counts: DataFrame containing the expression counts for the duplicated genes.
    :param rest_counts: DataFrame containing the expression counts for the rest of the genes.
    :return: Modified rest_counts DataFrame with merged expression counts for duplicated genes.
    """
    # Group rows in dupl_gene_counts by row names (after removing suffixes) and sum expression counts within each group
    dupl_gene_counts = dupl_gene_counts.groupby(dupl_gene_counts.index.str.replace(r'\\.[0-9]+$', '', regex=True)).sum()
    # Add dupl_gene_counts to rest_counts to combine expression counts for duplicated genes with rest of genes
    rest_counts = rest_counts.add(dupl_gene_counts, fill_value=0)
    # Remove suffixes from row names of rest_counts
    rest_counts.index = rest_counts.index.str.replace(r'\\.[0-9]+$', '', regex=True)
    return rest_counts


def Counts_filter(adata):
    # Filter out rows (genes) that contain lowercase letters in their names
    gene_name = adata.var_names
    adata = adata[np.array([not bool(re.search('[a-z]', x)) for x in gene_name]), :]
    gene_exp = np.sum(adata.X, axis=1)
    print(gene_exp)
    adata = adata[gene_exp > 0, :]
    # Filter out columns (cells) that have zero expression across all genes
    cell_exp = np.sum(adata.X, axis=0)
    print(cell_exp)
    adata = adata[:, cell_exp > 0]
    print(adata.var_names)
    return adata


def process_Xcounts(adata):
    """
    Processes an AnnData object by filtering out rows and columns and merging expression counts for duplicated genes.

    :param adata: AnnData object containing expression counts.
    :return: Modified AnnData object with filtered rows and columns and merged expression counts for duplicated genes.
    """
    # Filter out columns (cells) that contain lowercase letters in their names
    gene_name = adata.var_names
    adata = adata[:, np.array([not bool(re.search('[a-z]', x)) for x in gene_name])]
    
    # Filter out rows (genes) that have zero expression across all cells
    gene_exp = np.sum(adata.X, axis=1)
    adata = adata[gene_exp > 0, :]
    # Split expression counts into two DataFrames: one for duplicated genes and one for rest of genes
    fix_gene = adata.var_names.str.contains(r'\.[0-9]+$')
    # print(adata[fix_gene, :].var_names)
    gene_name = adata.var_names.str.replace(r'\.[0-9]+$', '', regex=True)
    dupl_gene_counts = adata[:, fix_gene].to_df()
    rest_counts = adata[:, ~fix_gene].to_df()
    # Merge expression counts for duplicated genes
    rest_counts = merge_duplicate_gene(dupl_gene_counts=dupl_gene_counts, rest_counts=rest_counts)
    
    return sc.AnnData(rest_counts)


"""\
Scanpy Merge
Merge tow or more accession samples via adata.concatenate()
Uses the implementation of *scikit-learn* [Pedregosa11]_.
.. version0.1
Parameters
----------
h5ad_path_list
    h5ad files path list.
merged_path
    output merged h5ad file path.
matrix_flag
    if separate R::Matrix files needed
    by default set to False.
dest
    separate R::Matrix files path
    Required if `matrix_flag=True` was passed.
Returns
-------

"""

def h5ad_merge(h5ad_files_path, merged_path, matrix_flag = False, dest = None):
    """
    Merges multiple h5ad files into one.

    :param h5ad_files_path: Path to a file containing a list of h5ad files to merge.
    :param merged_path: Path to save the merged h5ad file.
    :param matrix_flag: Whether to save the expression matrix as a .mtx file.
    :param dest: Destination directory to save the .mtx file.
    """
    adata_list = []
    h5ad_path_list = pd.read_csv(h5ad_files_path, sep=" ", header = None)[0].tolist()
    h5ad_path_list = [sub.replace('.rds', '.h5ad') for sub in h5ad_path_list]
    h5ad_path_list = [line for line in h5ad_path_list if re.search(r'BPCell', line) is None]
    total_cell = 0
    macro_group = 0
    # counter = 0
    for i in h5ad_path_list:
        adata = sc.read_h5ad(i)
        # new_adata = Counts_filter(adata)
        # # new_adata = (new_adata)
        # print(new_adata.var.index)
        adata_list.append(adata)
        total_cell += adata.n_obs
        # counter += 1
        if total_cell >= 900000:
            print(f"Current Num of Cells: {str(total_cell)}")
            adata = adata_list[0].concatenate(adata_list[1:], join='outer', batch_key='global')
            adata.obs['age'] = adata.obs['age'].astype(str)
            # View the structure and data types of the columns in adata.obs
            print(adata.obs.info())
            print(adata.var.index)
            group_path = str(merged_path.replace('.h5ad', '_MacroGroup'+str(macro_group)+'.h5ad'))
            
            print(group_path)
            # store adata into a h5ad format
            adata.write(group_path)
            meta_path = re.sub(".h5ad", "_meta.csv", group_path)
            adata.obs.to_csv(meta_path, sep = ",")
            print(f"Successfully merge MacroGroup{str(macro_group)}")
            total_cell = 0
            macro_group += 1
            adata_list.clear()
            

    print("Merge the rest h5ad files")
    print(f"Current Num of Cells: {str(total_cell)}")
    # Integrate the rest h5ad object
    adata = adata_list[0].concatenate(adata_list[1:], join='outer', batch_key='global')
    adata.obs['age'] = adata.obs['age'].astype(str)
    print(adata.obs.info())

    group_path = str(merged_path.replace('.h5ad', '_MacroGroup'+str(macro_group)+'.h5ad'))
    print(group_path)
    # store adata into a h5ad format
    adata.write(group_path)
    meta_path = re.sub(".h5ad", "_meta.csv", group_path)
    adata.obs.to_csv(meta_path, sep = ",")


    if matrix_flag == True:
        pd.DataFrame(adata.var.index).to_csv(os.path.join(dest, "genes.tsv" ),   sep = "\t")
        pd.DataFrame(adata.obs.index).to_csv(os.path.join(dest, "barcodes.tsv"), sep = "\t")
        scipy.io.mmwrite(os.path.join(dest, "matrix.mtx"), adata.X)
    print("Successfully Merge h5ad file lists\n")



if __name__ == "__main__":
    # print(f"Arguments count: {len(sys.argv)}")
    h5ad_merge(h5ad_files_path = sys.argv[1], merged_path = sys.argv[2], matrix_flag = sys.argv[3])
    try:
        print("Hello")
        

    except:
        print("Error with h5ad merge\n")