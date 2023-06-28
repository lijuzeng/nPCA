import numpy as np
import scanpy as sc
import pandas as pd


def rawdata2data(file):
    # rawdata = sc.read_10x_mtx(file)                                                     # read-in .mtx data
    rawdata = sc.read_text(file, delimiter='\t', first_column_names=True)                 # read-in .csv data
    rawdata.var_names_make_unique()
    print(rawdata)

    sc.pp.filter_cells(rawdata, min_genes=200)                  # Remove cells expressing less than 200 genes
    sc.pp.filter_genes(rawdata, min_cells=3)                    # Remove genes expressing below 3 cells
    print(rawdata)

    # Calculate mitochondrial gene expression level/all gene expression level for each cell
    # mito_gene = rawdata.var_names.str.startswith('MT-')       # human
    mito_gene = rawdata.var_names.str.startswith('Mt-')         # mouse
    rawdata.obs['percent_mito'] = \
        np.ravel(np.mat(rawdata[:, mito_gene].X.sum(axis=1))) / np.ravel(np.mat(rawdata.X.sum(axis=1)))
    rawdata = rawdata[rawdata.obs.n_genes < 2500, :]                # Obtain cell samples with less than 2500 expressed genes
    rawdata = rawdata[rawdata.obs.percent_mito < 0.05, :]           # Obtain cell samples with less than 5% mitochondrial genes
    print(rawdata)

    sc.pp.normalize_total(rawdata, target_sum=1e4)
    sc.pp.log1p(rawdata)
    sc.pp.highly_variable_genes(rawdata, n_top_genes=1000)
    # sc.pp.highly_variable_genes(rawdata, min_mean=0.05, max_mean=0.2)
    print(rawdata)
    rawdata = rawdata[:, rawdata.var["highly_variable"]]
    print(rawdata)
    rawdata.to_df().to_csv("./Baron_after_filter.csv", sep='\t')


if __name__ == "__main__":
    rawdata2data("/home/ljz/project/npca/data/single_cell/Baron/Baron_mouse2/Baron_mouse2.data")
    # rawdata2data("/home/ljz/project/npca/data/single_cell/3k_Human_Squamous_Cell_Lung_Carcinoma_DTCs/filtered_feature_bc_matrix")
