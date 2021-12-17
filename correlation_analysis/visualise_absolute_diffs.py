from urllib.error import HTTPError
from warnings import resetwarnings
import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Activate rpy2 
pandas2ri.activate()
readRDS = robjects.r['readRDS']

# Read UID: Gene symbol mapper

def get_dict(gene_codes):
    df = readRDS(gene_codes)
    genedict={df.iloc[i,0]:df.iloc[i,1] for i in range(len(df))}
    genedict["ENSG00000133895"]="MEN1"
    return genedict

# Import and clean correlation matrices

def get_matrices(cancer, healthy):
    df_c = readRDS(cancer)
    df_h= readRDS(healthy)
    index = pd.Series([i[:-3] for i in df_h.index.to_list()])
    df_h=df_h.set_index(index)
    df_h=df_h.sort_index()
    df_c=df_c.sort_index()
    new_index = [i for i in df_c.index if i in df_h.index]
    df_c=df_c[df_c.index.isin(new_index)]
    return df_c, df_h

# Merge datasets and obtain absolute differences

def merge_sort(df_c, df_h):
    merged=df_h.join(df_c)
    merged.columns=['Adrenal Gland_h','Bladder_h','Brain_h','Lung_h','Thyroid_h',
    'Adrenal Gland_c','Bladder_c','Brain_c','Lung_c','Thyroid_c']
    idx = merged.index.tolist()
    out=[[],[],[],[],[],[]]
    for i in range(len(merged)):
        out[0].append(merged.index.tolist()[i])
        out[1].append(np.abs(merged.iloc[i,0]-merged.iloc[i,5]))
        out[2].append(np.abs(merged.iloc[i,1]-merged.iloc[i,6]))
        out[3].append(np.abs(merged.iloc[i,2]-merged.iloc[i,7]))
        out[4].append(np.abs(merged.iloc[i,3]-merged.iloc[i,8]))
        out[5].append(np.abs(merged.iloc[i,4]-merged.iloc[i,9]))
    diff = pd.DataFrame(np.asarray(out[1:])).T
    diff.columns=['Adrenal Gland','Bladder','Brain','Lung','Thyroid']
    diff.index=idx
    diff['sum']=diff.sum(axis=1)
    diff=diff.dropna(axis=0)
    diff=diff.sort_values('sum',ascending=False)
    top = diff.iloc[:5,:]
    bottom = diff.iloc[-5:,:]
    merged = pd.concat([top,bottom])
    m=merged.to_numpy()

    return merged, m


def plot_diffs(merged, m, genedict):
    codes=[genedict[x] for x in merged.index.tolist()]

    fig, ax = plt.subplots(figsize=(5,5.5))
    ax.grid(False)
    im=ax.imshow(m)
    xlabs=merged.columns.tolist()
    ylabs=codes
    ax.set_xticks(np.arange(len(xlabs)))
    ax.set_xticklabels(xlabs)
    ax.set_yticks(np.arange(len(ylabs)))
    ax.set_yticklabels(ylabs)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")
    plt.setp(ax.get_yticklabels(), rotation=0, ha="right",
            rotation_mode="anchor")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    plt.show()

def run(gene_codes,cancer_rds,healthy_rds):
    genedict = get_dict(gene_codes)
    df_c, df_h = get_matrices(cancer_rds, healthy_rds)
    merged, m= merge_sort(df_c,df_h)
    plot_diffs(merged,m,genedict)

if __name__== "__main__":
    gene_codes='DDR_epi_genes_names.rds'
    cancer = "correlation_analysis/unfiltered_cancer_correlation_matrix.rds"
    healthy = "correlation_analysis/unfiltered_gtex_correlation_matrix.rds"
    run(gene_codes,cancer,healthy)

