import numpy as np
import scanpy as sc
import pandas as pd


def rawdata2data(file):
    # rawdata = sc.read_10x_mtx(file)                             # 读入 .mtx 数据
    rawdata = sc.read_text(file, delimiter='\t', first_column_names=True)                 # 读入 .csv 数据
    rawdata.var_names_make_unique()
    print(rawdata)

    sc.pp.filter_cells(rawdata, min_genes=200)                  # 去除表达基因200以下的细胞
    sc.pp.filter_genes(rawdata, min_cells=3)                    # 去除在3个细胞以下表达的基因
    print(rawdata)

    # 对每一个细胞计算 线粒体基因表达量 / 所有基因表达量
    # mito_gene = rawdata.var_names.str.startswith('MT-')         # human
    mito_gene = rawdata.var_names.str.startswith('Mt-')         # mouse
    rawdata.obs['percent_mito'] = \
        np.ravel(np.mat(rawdata[:, mito_gene].X.sum(axis=1))) / np.ravel(np.mat(rawdata.X.sum(axis=1)))
    rawdata = rawdata[rawdata.obs.n_genes < 2500, :]                # 获取表达基因数在 2500 以下的细胞样本
    rawdata = rawdata[rawdata.obs.percent_mito < 0.05, :]           # 获取线粒体基因占比在 5% 以下的细胞样本
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
