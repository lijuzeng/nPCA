import numpy as np
import matplotlib.gridspec as gridspec
from sklearn.decomposition import FastICA as ICA
from sklearn.manifold import TSNE
from sklearn.manifold import MDS
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.metrics import adjusted_rand_score
import matplotlib.pyplot as plt

file_path = '/home/ljz/project/npca/result/single_cell/Baron/Baron_mouse2/1000/'
raw_data = np.loadtxt("/home/ljz/project/npca/data/single_cell/Baron/Baron_mouse2/1000/Baron_mouse2_1000.data", delimiter='\t', encoding='utf-8')
pca_data = np.loadtxt(file_path+'pc.pca', delimiter='\t', encoding='utf-8')
npca_data = np.loadtxt(file_path+'pc.npca', delimiter='\t', encoding='utf-8')
true_label = np.loadtxt("/home/ljz/project/npca/data/single_cell/Baron/Baron_mouse2/Baron_mouse2_label.data")

# run ICA
ica = ICA(n_components=2, random_state=0, whiten='unit-variance').fit_transform(raw_data)

# run MDS
mds = MDS(n_components=2, random_state=0).fit_transform(raw_data)

def true(s, file):
    # data of  ica \ mds \ pca \ npca 
    ix = ica[:, 0]
    mx = mds[:, 0]
    x = pca_data[:, 0]
    nx = npca_data[:, 0]
    iy = ica[:, 1]
    my = mds[:, 1]
    y = pca_data[:, 1]
    ny = npca_data[:, 1]
    # construct the figure
    fig = plt.figure(figsize=(14, 10))
    grid = gridspec.GridSpec(2, 8)
    ax1 = fig.add_subplot(grid[0, 1:4])
    ax2 = fig.add_subplot(grid[0, 5:8])
    ax3 = fig.add_subplot(grid[1, 1:4])
    ax4 = fig.add_subplot(grid[1, 5:8])
    # set the title
    ax1.set_title('ica')
    ax2.set_title('mds')
    ax3.set_title('pca')
    ax4.set_title('npca')
    # set abscissa
    ax1.set_xlabel('ica_1')
    ax2.set_xlabel('mds_1')
    ax3.set_xlabel('pc_1')
    ax4.set_xlabel('pc_1')
    # set ordinate
    ax1.set_ylabel('ica_2')
    ax2.set_ylabel('mds_2')
    ax3.set_ylabel('pc_2')
    ax4.set_ylabel('pc_2')
    # draw
    plot1 = ax1.scatter(ix, iy, s=s, c=true_label, cmap='tab20')
    plot2 = ax2.scatter(mx, my, s=s, c=true_label, cmap='tab20')
    plot3 = ax3.scatter(x, y, s=s, c=true_label, cmap='tab20')
    plot4 = ax4.scatter(nx, ny, s=s, c=true_label, cmap='tab20')
    # add legend
    a, b = plot1.legend_elements()
    b = ['$\\mathdefault{activated\_stellate}$', 
        '$\\mathdefault{alpha}$', 
        '$\\mathdefault{B\_cell}$', 
        '$\\mathdefault{beta}$', 
        '$\\mathdefault{delta}$', 
        '$\\mathdefault{ductal}$', 
        '$\\mathdefault{endothelial}$',  
        '$\\mathdefault{gamma}$', 
        '$\\mathdefault{immune\_other}$', 
        '$\\mathdefault{macrophage}$', 
        '$\\mathdefault{quiescent\_stellate}$', 
        '$\\mathdefault{T\_cell}$']
    ax1.legend(a, b, loc=3, bbox_to_anchor=(-0.7, -0.5), borderaxespad=0)
    # legend2 = ax2.legend(*plot1.legend_elements())
    # legend2 = ax2.legend(*plot2.legend_elements())
    # legend3 = ax3.legend(*plot3.legend_elements())
    # legend4 = ax4.legend(*plot4.legend_elements())
    plt.savefig(file+'true_label.png', dpi=300)
    plt.show()

def kmean_2d(k, s, file):
    km_ica = KMeans(n_clusters=k, random_state=0).fit(ica)
    kmica_label = km_ica.labels_
    km_mds = KMeans(n_clusters=k, random_state=0).fit(mds)
    kmmds_label = km_mds.labels_
    km_pca = KMeans(n_clusters=k, random_state=0).fit(pca_data)
    kmpca_label = km_pca.labels_
    km_npca = KMeans(n_clusters=k, random_state=0).fit(npca_data)
    kmnpca_label = km_npca.labels_

    # data of  ica \ mds \ pca \ npca 
    ix = ica[:, 0]
    mx = mds[:, 0]
    x = pca_data[:, 0]
    nx = npca_data[:, 0]
    iy = ica[:, 1]
    my = mds[:, 1]
    y = pca_data[:, 1]
    ny = npca_data[:, 1]
    # construct the figure
    fig = plt.figure(figsize=(14, 10))
    grid = gridspec.GridSpec(2, 8)
    ax1 = fig.add_subplot(grid[0, 1:4])
    ax2 = fig.add_subplot(grid[0, 5:8])
    ax3 = fig.add_subplot(grid[1, 1:4])
    ax4 = fig.add_subplot(grid[1, 5:8])
    # set the title
    ax1.set_title('kmean of ica')
    ax2.set_title('kmean of mds')
    ax3.set_title('kmean of pca')
    ax4.set_title('kmean of npca')
    # set abscissa
    ax1.set_xlabel('ica_1')
    ax2.set_xlabel('mds_1')
    ax3.set_xlabel('pc_1')
    ax4.set_xlabel('pc_1')
    # set ordinate
    ax1.set_ylabel('ica_2')
    ax2.set_ylabel('mds_2')
    ax3.set_ylabel('pc_2')
    ax4.set_ylabel('pc_2')
    # draw
    plot1 = ax1.scatter(ix, iy, s=s, c=kmica_label, cmap='tab20')
    plot2 = ax2.scatter(mx, my, s=s, c=kmmds_label, cmap='tab20')
    plot3 = ax3.scatter(x, y, s=s, c=kmpca_label, cmap='tab20')
    plot4 = ax4.scatter(nx, ny, s=s, c=kmnpca_label, cmap='tab20')
    # add legend
    a, b = plot1.legend_elements()
    b = ['$\\mathdefault{activated\_stellate}$', 
        '$\\mathdefault{alpha}$', 
        '$\\mathdefault{B\_cell}$', 
        '$\\mathdefault{beta}$', 
        '$\\mathdefault{delta}$', 
        '$\\mathdefault{ductal}$', 
        '$\\mathdefault{endothelial}$',  
        '$\\mathdefault{gamma}$', 
        '$\\mathdefault{immune\_other}$', 
        '$\\mathdefault{macrophage}$', 
        '$\\mathdefault{quiescent\_stellate}$', 
        '$\\mathdefault{T\_cell}$']
    ax1.legend(a, b, loc=3, bbox_to_anchor=(-0.7, -0.5), borderaxespad=0)
    # legend2 = ax2.legend(*plot2.legend_elements())
    # legend3 = ax3.legend(*plot3.legend_elements())
    # legend4 = ax4.legend(*plot4.legend_elements())
    plt.savefig(file+'kmean.png', dpi=300)
    plt.show()
    # calculate silhouette_score
    ica_ss = silhouette_score(ica, kmica_label)
    mds_ss = silhouette_score(mds, kmmds_label)
    pca_ss = silhouette_score(pca_data, kmpca_label)
    npca_ss = silhouette_score(npca_data, kmnpca_label)
    print('--------------------------------------')
    print("silhouette_score of ica_kmean = ", ica_ss)
    print("silhouette_score of mds_kmean = ", mds_ss)
    print("silhouette_score of pca_kmean = " , pca_ss)
    print("silhouette_score of npca_kmean = ", npca_ss)
    print('--------------------------------------')
    # calculate adjusted_rank_socre
    ica_ari = adjusted_rand_score(true_label, kmica_label)
    mds_ari = adjusted_rand_score(true_label, kmmds_label)
    pca_ari = adjusted_rand_score(true_label, kmpca_label)
    npca_ari = adjusted_rand_score(true_label, kmnpca_label)
    print('--------------------------------------')
    print("adjusted_rand_score of ica_kmean = ", ica_ari)
    print("adjusted_rand_score of mds_kmean = ", mds_ari)
    print("adjusted_rand_score of pca_kmean = " , pca_ari)
    print("adjusted_rand_score of npca_kmean = ", npca_ari)
    print('--------------------------------------')

def tsne_2d(s, file):
    tsne_ica = TSNE(n_components=2, random_state=0, learning_rate='auto', init='random').fit(ica)
    tsneica_embedding = tsne_ica.embedding_
    tsne_mds = TSNE(n_components=2, random_state=0, learning_rate='auto', init='random').fit(mds)
    tsnemds_embedding = tsne_mds.embedding_
    tsne_pca = TSNE(n_components=2, random_state=0, learning_rate='auto', init='random').fit(pca_data)
    tsnepca_embedding = tsne_pca.embedding_
    tsne_npca = TSNE(n_components=2, random_state=0, learning_rate='auto', init='random').fit(npca_data)
    tsnenpca_embedding = tsne_npca.embedding_

    # data
    ix= tsneica_embedding[:, 0]
    mx = tsnemds_embedding[:, 0]
    x = tsnepca_embedding[: ,0]
    nx = tsnenpca_embedding[:, 0]
    iy = tsneica_embedding[:, 1]
    my = tsnemds_embedding[:, 1]
    y = tsnepca_embedding[:, 1]
    ny = tsnenpca_embedding[:, 1]
    # construct the figure
    fig = plt.figure(figsize=(14, 10))
    grid = gridspec.GridSpec(2, 8)
    ax1 = fig.add_subplot(grid[0, 1:4])
    ax2 = fig.add_subplot(grid[0, 5:8])
    ax3 = fig.add_subplot(grid[1, 1:4])
    ax4 = fig.add_subplot(grid[1, 5:8])
    # set the title
    ax1.set_title('t-SNE of ica')
    ax2.set_title('t-SNE of mds')
    ax3.set_title('t-SNE of pca')
    ax4.set_title('t-SNE of npca')
    # set abscissa 
    ax1.set_xlabel('tSNE_1')
    ax2.set_xlabel('tSNE_1')
    ax3.set_xlabel('tSNE_1')
    ax4.set_xlabel('tSNE_1')
    # set ordinate
    ax1.set_ylabel('tSNE_2')
    ax2.set_ylabel('tSNE_2')
    ax3.set_ylabel('tSNE_2')
    ax4.set_ylabel('tSNE_2')
    # draw
    plot1 = ax1.scatter(ix, iy, s=s, c=true_label, cmap='tab20')
    plot2 = ax2.scatter(mx, my, s=s, c=true_label, cmap='tab20')
    plot3 = ax3.scatter(x, y, s=s, c=true_label, cmap='tab20')
    plot4 = ax4.scatter(nx, ny, s=s, c=true_label, cmap='tab20')
    # add legend
    a, b = plot1.legend_elements()
    b = ['$\\mathdefault{activated\_stellate}$', 
        '$\\mathdefault{alpha}$', 
        '$\\mathdefault{B\_cell}$', 
        '$\\mathdefault{beta}$', 
        '$\\mathdefault{delta}$', 
        '$\\mathdefault{ductal}$', 
        '$\\mathdefault{endothelial}$',  
        '$\\mathdefault{gamma}$', 
        '$\\mathdefault{immune\_other}$', 
        '$\\mathdefault{macrophage}$', 
        '$\\mathdefault{quiescent\_stellate}$', 
        '$\\mathdefault{T\_cell}$']
    ax1.legend(a, b, loc=3, bbox_to_anchor=(-0.7, -0.5), borderaxespad=0)
    # legend2 = ax2.legend(*plot2.legend_elements())
    # legend3 = ax3.legend(*plot3.legend_elements())
    # legend4 = ax4.legend(*plot4.legend_elements())
    plt.savefig(file+'tsne.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    true(1, file_path)
    kmean_2d(12, 1, file_path)
    tsne_2d(1, file_path)
