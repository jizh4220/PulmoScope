# %%
import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
# import argparse

# %%

# %%


def refine_metanames(meta):
    meta.rename(columns={
        "orig.ident": "accession",
        "gse_alias": "experiment",
        "DISEASE": "disease",
        "TISSUE": "tissue",
        "celltype": "celltype"
        # "macro_celltype": "celltype"
    }, inplace=True)
    return meta

def clean_metanoise(meta):
    meta = refine_metanames(meta)
    meta_param = ["disease", "tissue", "experiment", "celltype", "groups", "accession"]
    meta = meta[meta_param]
    meta.replace(["", "NA", None], "unknown", inplace=True)
    return meta


def fraction_param_table(adata, threshold=0.5, param='param', goi='SCGB1A1'):
    # preprocess meta data names
    adata.obs['param'] = adata.obs[param]
    meta = adata.obs
    meta_param = ["disease", "tissue", "experiment", "celltype", "groups", "accession", 'param']
    frac_meta = meta.loc[:, meta.columns.intersection(meta_param)].drop_duplicates(subset=['accession', param])
    frac_meta['idx'] = frac_meta[param].astype(str) + frac_meta['accession'].astype(str)

    df = adata[:, goi].to_df()
    df['accession'] = adata.obs['accession']
    df['param'] = adata.obs['param']
    df['disease'] = adata.obs['disease']
    # Create a DataFrame for 'acc_table'
    acc_table = df.groupby(['accession', 'param']).size().reset_index(name='Freq')
    acc_table['idx'] = acc_table['param'].astype(str) + acc_table['accession'].astype(str)
    acc_table = acc_table[acc_table['Freq'] > 0]

    frac_meta = pd.DataFrame()
    frac_meta['idx'] = df['param'].astype(str) + df['accession'].astype(str)
    frac_meta['disease'] = df['disease']
    # Merge frac_meta and acc_table on 'idx'
    frac_meta = pd.merge(frac_meta, acc_table, how='left', left_on='idx', right_on='idx')
    # Fill NaN values with 0
    frac_meta['Freq'] = frac_meta['Freq'].fillna(0)
    # Rename 'Freq' to 'cell_counts'
    frac_meta.rename(columns={'Freq': 'cell_counts'}, inplace=True)
    cluster_frac = []
    # Every gene expressed counts
    for i in range(len(goi)):
        # current gene cluster-wise expression
        tmp = df.iloc[:, i]
        exp_cells = df[tmp >= threshold]
        
        # Create a DataFrame for 'exp_table'
        exp_table = exp_cells.groupby(['accession', 'param']).size().reset_index(name='Freq')
        exp_table = exp_table[exp_table['Freq'] > 0]
        goi_frac = frac_meta.copy()
        goi_frac['gene'] = goi[i]
        goi_frac['expressed_counts'] = 0
        goi_frac[param] = goi_frac['param']
        goi_frac.loc[exp_table.index, 'expressed_counts'] = exp_table['Freq'].values
        goi_frac['expressed_ratio'] = goi_frac['expressed_counts'] / goi_frac['cell_counts']
        cluster_frac.append(goi_frac)

    frac_table = pd.concat(cluster_frac)
    frac_table.drop(columns=['idx'], inplace=True)
    return frac_table

def fraction_plot(frac_table, goi='SCGB1A1', threshold=0.5, param='disease', floor=50):
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Assuming 'frac_table' is your DataFrame and 'goi' is your list of genes of interest
    frac_table['param'] = frac_table[param]
    sample_names = frac_table['param'].unique()
    ps = []

    for gene in goi:
        cur_frac = frac_table[frac_table['gene'] == gene]
        all_sample_table = [len(cur_frac['accession'][cur_frac['param'] == sample].unique()) for sample in sample_names]
        cur_frac = cur_frac[(cur_frac['cell_counts'] >= floor) & (cur_frac['expressed_ratio'] > 0)]
        
        if cur_frac.empty:
            print(f"Current gene: {gene} is not found in Fraction Table")
            continue
        
        cur_frac['param_ct'] = cur_frac['celltype'].astype(str) + "_" + cur_frac['param'].astype(str)
        sample_table = [len(cur_frac['accession'][cur_frac['param'] == sample].unique()) for sample in sample_names]
        
        plt.figure(figsize=(10, 6))
        sns.boxplot(x='param_ct', y='expressed_ratio', data=cur_frac, hue='param')
        sns.swarmplot(x='param_ct', y='expressed_ratio', data=cur_frac, color=".25")
        
        title = f"Fraction plot of GOI: {gene}"
        subtitle = (
            f"Total {sample_names[0]} Sample Number: {all_sample_table[0]}\n"
            f"Total {sample_names[0]} Sample Number Pass Cutoff: {sample_table[0]}\n"
            f"Total {sample_names[1]} Sample Number: {all_sample_table[1]}\n"
            f"Total {sample_names[1]} Sample Number Pass Cutoff: {sample_table[1]}\n"
            f"Total Cell Counts: {frac_table['cell_counts'].nunique()}\n"
            f"Pass Cutoff Cell Counts: {round(sum(cur_frac['cell_counts'] * cur_frac['expressed_ratio']))}\n"
            f"Expression Cutoff: {threshold}\n"
            f"Minimum Cell Number: {floor}"
        )
        
        plt.title(title, y=1.5, loc='left')  # Align the title to the left
        plt.text(x=0, y=1.5, s=subtitle, ha='left', va='top', transform=plt.gca().transAxes)  # Add the subtitle using plt.text()
        plt.xlabel('Genes of Interest')
        plt.ylabel('Fraction of Cells Expressing Target Genes')
        plt.xticks(rotation=45)
        
        plt.ylim(0, cur_frac['expressed_ratio'].mean() + 1)  # Adjust 0.1 to provide some padding
        # Append the current figure to the list

        # Adjust the top margin
        plt.subplots_adjust(top=0.6)
        plt.savefig(f'{gene}_FractionPlot.png', bbox_inches='tight')
        ps.append(plt.gcf())
    return ps



# %%
# clean_meta = refine_metanames(adata.obs)
# adata.obs = clean_meta
# frac_table = fraction_param_table(adata, threshold=0.5, param='celltype', goi=['Ddr2'])

# # %%
# ps = fraction_plot(frac_table, goi=['Ddr2'], threshold=0.5, param='disease', floor=50)
# # Display all plots
# for fig in ps:
#     fig.show()

