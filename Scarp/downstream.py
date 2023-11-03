import sys
sys.path.append('./Scarp/')

from SCARP_help_func import *
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
from sklearn.metrics import confusion_matrix
import seaborn as sns
import scipy
from sklearn import metrics


'''
#########################################################################################################
#                                        downstream analysis
#########################################################################################################
'''


def construct_adata(data_name, mat, Cells, labels,
                    cluster_num, random_state,
                    plot_fig=False, plot_file=None,
                    plot_confusion=False, conf_file=None):
    Cells_df = pd.DataFrame(mat,
                            index=Cells,
                            columns=['Peaks' + str(i + 1) for i in range(mat.shape[1])])

    adata_cell = sc.AnnData(Cells_df)
    adata_cell.var_names_make_unique()
    adata_cell.obs['label'] = labels
    sc.pp.neighbors(adata_cell, use_rep='X', random_state=random_state)

    [warning, adata_cell] = getNClusters_Louvain(adata_cell, cluster_num, range_max=3)
    if warning == 0:
        print('Wrong for Louvain! ')

    [warning, adata_cell] = getNClusters_Leiden(adata_cell, cluster_num, range_max=3)
    if warning == 0:
        print('Wrong for Leiden! ')

    if data_name == 'GM12878vsHEK' or data_name == 'GM12878vsHL':
        apart_mixed_label = adata_cell.obs['label'].iloc[np.where(adata_cell.obs['label'] != 'Mixed')].astype(
            'str').astype('category')
        apart_mixed_louvain = adata_cell.obs['louvain'].iloc[np.where(adata_cell.obs['label'] != 'Mixed')].astype(
            'str').astype('category')
        apart_mixed_leiden = adata_cell.obs['leiden'].iloc[np.where(adata_cell.obs['label'] != 'Mixed')].astype(
            'str').astype('category')
    else:
        apart_mixed_label = adata_cell.obs['label']
        apart_mixed_louvain = adata_cell.obs['louvain']
        apart_mixed_leiden = adata_cell.obs['leiden']

    score = pd.DataFrame({'ARI': [adjusted_rand_score(apart_mixed_louvain, apart_mixed_label),
                                  adjusted_rand_score(apart_mixed_leiden, apart_mixed_label)],
                          'NMI': [normalized_mutual_info_score(apart_mixed_louvain, apart_mixed_label),
                                  normalized_mutual_info_score(apart_mixed_leiden, apart_mixed_label)]},
                         index=['louvain', 'leiden'])

    if plot_fig:
        sc.tl.umap(adata_cell, random_state=random_state)
        fig1, ax = plt.subplots(figsize=(3, 3))
        sc.pl.umap(adata_cell, color='label', title='labeled', s=20, ax=ax,
                   save='_' + data_name + '_' + plot_file + '_labeled.svg')
        sc.tl.tsne(adata_cell, use_rep='X', random_state=random_state)
        fig2, ax = plt.subplots(figsize=(3, 3))
        sc.pl.tsne(adata_cell, color='label', title='labeled', s=20, ax=ax,
                   save='_' + data_name + '_' + plot_file + '_labeled.svg')
        fig3, ax = plt.subplots(figsize=(3, 3))
        sc.pl.umap(adata_cell, color='louvain', title='Louvain', s=20, ax=ax,
                   save='_' + data_name + '_' + plot_file + '_louvain.svg')
        fig4, ax = plt.subplots(figsize=(3, 3))
        sc.pl.umap(adata_cell, color='leiden', title='Leiden', s=20, ax=ax,
                   save='_' + data_name + '_' + plot_file + '_leiden.svg')

    if plot_confusion:
        if data_name == 'GM12878vsHL':
            name = dict(zip(['GM12878', 'HL60', 'Mixed'], range(3)))
        elif data_name == 'GM12878vsHEK':
            name = dict(zip(['GM12878', 'HEK293T', 'Mixed'], range(3)))
        else:
            name = dict(zip(set(labels), range(cluster_num)))

        Cluster_louvain_mat = confusion_matrix(adata_cell.obs['louvain'].astype(int), [name[i] for i in labels])
        Cluster_leiden_mat = confusion_matrix(adata_cell.obs['leiden'].astype(int), [name[i] for i in labels])

        if data_name == 'GM12878vsHEK' or data_name == 'GM12878vsHL':
            Cluster_louvain_mat = pd.DataFrame(Cluster_louvain_mat[:-1, :],
                                               columns=list(name.keys()),
                                               index=range(cluster_num))
            Cluster_leiden_mat = pd.DataFrame(Cluster_leiden_mat[:-1, :],
                                              columns=list(name.keys()),
                                              index=range(cluster_num))
        else:
            Cluster_louvain_mat = pd.DataFrame(Cluster_louvain_mat, columns=list(name.keys()),
                                               index=range(cluster_num))
            Cluster_leiden_mat = pd.DataFrame(Cluster_leiden_mat, columns=list(name.keys()),
                                              index=range(cluster_num))

        Cluster_louvain_mat = confusion_sort(Cluster_louvain_mat)
        Cluster_leiden_mat = confusion_sort(Cluster_leiden_mat)
        Cluster_louvain_mat.to_csv(conf_file + data_name + '_louvain_confusion_mat.csv')
        Cluster_leiden_mat.to_csv(conf_file + data_name + '_leiden_confusion_mat.csv')

        plt.figure(figsize=(4, 4))
        sns.heatmap(Cluster_louvain_mat, annot=True, cmap="Blues", fmt='.20g', cbar=False)
        plt.xticks(rotation=45, horizontalalignment="center")
        plt.ylabel('Predicted')
        plt.xlabel('True')
        plt.title(data_name + '(louvain)', size=10)
        plt.savefig(conf_file + data_name + '_louvain_confusion_mat.svg', bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(4, 4))
        sns.heatmap(Cluster_leiden_mat, annot=True, cmap="Blues", fmt='.20g', cbar=False)
        plt.xticks(rotation=45, horizontalalignment="center")
        plt.ylabel('Predicted')
        plt.xlabel('True')
        plt.title(data_name + '(leiden)', size=10)
        plt.savefig(conf_file + data_name + '_leiden_confusion_mat.svg', bbox_inches='tight')
        plt.close()

    return score


def compute_score(data_name, diffusion_mat, cluster_num, labels, Cells,
                  kept_comp=10, random_state=1,
                  plot_fig=False, plot_file=None,
                  plot_confusion=False, conf_file=None):
    cell_mat = SCARP_cell_embedding(diffusion_mat, kept_comp)
    score = construct_adata(data_name=data_name,
                            mat=cell_mat,
                            Cells=Cells,
                            labels=labels,
                            cluster_num=cluster_num,
                            random_state=random_state,
                            plot_fig=plot_fig,
                            plot_file=plot_file,
                            plot_confusion=plot_confusion,
                            conf_file=conf_file)
    return score


def confusion_sort(confusion_mat):
    cluster_num = confusion_mat.shape[0]
    mat_sort = np.sort(np.array(confusion_mat).flatten())
    fix = []
    for i in range(1, confusion_mat.shape[0]*confusion_mat.shape[1]):
        temp_max = mat_sort[-i]
        pre, tru = np.where(confusion_mat == temp_max)

        for t in range(tru.shape[0]):
            if (tru[t] not in fix) & (pre[t] not in fix):
                temp = confusion_mat.iloc[tru[t], :].copy()
                confusion_mat.iloc[tru[t], :] = confusion_mat.iloc[pre[t], :].copy()
                confusion_mat.iloc[pre[t], :] = temp
                fix.append(tru[t])

        if len(fix) == cluster_num:
            break

    return confusion_mat


def getNClusters_Louvain(adata, n_cluster, range_min=0, range_max=3, max_steps=50, random_seed=1):
    np.random.seed(random_seed)
    temp_step = 0
    temp_min = float(range_min)
    temp_max = float(range_max)
    while temp_step < max_steps:
        temp_resolution = temp_min + ((temp_max - temp_min) / 2)
        sc.tl.louvain(adata, resolution=temp_resolution)
        temp_clusters = adata.obs['louvain'].nunique()

        if temp_clusters > n_cluster:
            temp_max = temp_resolution
        elif temp_clusters < n_cluster:
            temp_min = temp_resolution
        else:
            return [1, adata]
        temp_step += 1

    print('Cannot find the number of clusters')
    print('Clustering solution from last iteration is used:' +
          str(temp_clusters) + ' at resolution ' + str(temp_resolution))
    sc.tl.louvain(adata, resolution=temp_resolution)
    adata.obs['louvain'].nunique()
    return [0, adata]


def getNClusters_Leiden(adata, n_cluster, range_min=0, range_max=3, max_steps=50, random_seed=1):
    np.random.seed(random_seed)
    temp_step = 0
    temp_min = float(range_min)
    temp_max = float(range_max)
    while temp_step < max_steps:
        temp_resolution = temp_min + ((temp_max - temp_min) / 2)
        sc.tl.leiden(adata, resolution=temp_resolution)
        temp_clusters = adata.obs['leiden'].nunique()

        if temp_clusters > n_cluster:
            temp_max = temp_resolution
        elif temp_clusters < n_cluster:
            temp_min = temp_resolution
        else:
            return [1, adata]
        temp_step += 1

    print('Cannot find the number of clusters')
    print('Clustering solution from last iteration is used:' +
          str(temp_clusters) + ' at resolution ' + str(temp_resolution))
    sc.tl.leiden(adata, resolution=temp_resolution)
    adata.obs['leiden'].nunique()
    return [0, adata]


def visualization(data_name, cell_embedding, cluster_num, labels, Cells,
                  random_state=1, save_fig=False, plot_confusion=False,
                  save_file=None):
    Cells_df = pd.DataFrame(cell_embedding,
                            index=Cells,
                            columns=['Peaks' + str(i + 1) for i in range(cell_embedding.shape[1])])

    adata_cell = sc.AnnData(Cells_df)
    adata_cell.var_names_make_unique()
    adata_cell.obs['label'] = labels

    # ==========================================================================================
    #                                          plot
    # ==========================================================================================
    sns.set_style("white")
    sc.pp.neighbors(adata_cell, use_rep='X', random_state=random_state)
    sc.tl.umap(adata_cell, random_state=random_state)
    sc.tl.tsne(adata_cell, use_rep='X', random_state=random_state)
    # labeled results
    fig1, ax = plt.subplots(figsize=(5, 4))
    sc.pl.umap(adata_cell, color='label', title='labeled', s=20, ax=ax, color_map='Blues')
    plt.subplots_adjust(left=0.1, right=0.7, top=0.9, bottom=0.1)
    if save_fig:
        plt.savefig(save_file + 'umap.csv')

    fig2, ax = plt.subplots(figsize=(5, 4))
    sc.pl.tsne(adata_cell, color='label', title='labeled', s=20, ax=ax)
    plt.subplots_adjust(left=0.1, right=0.7, top=0.9, bottom=0.1)
    if save_fig:
        plt.savefig(save_file + 'tsne.csv')

    [warning, adata_cell] = getNClusters_Louvain(adata_cell, cluster_num, range_max=3)
    if warning == 0:
        print('Wrong for Louvain! ')

    [warning, adata_cell] = getNClusters_Leiden(adata_cell, cluster_num, range_max=3)
    if warning == 0:
        print('Wrong for Leiden! ')

    if data_name == 'GM12878vsHEK' or data_name == 'GM12878vsHL':
        apart_mixed_label = adata_cell.obs['label'].iloc[np.where(adata_cell.obs['label'] != 'Mixed')].astype(
            'str').astype('category')
        apart_mixed_louvain = adata_cell.obs['louvain'].iloc[np.where(adata_cell.obs['label'] != 'Mixed')].astype(
            'str').astype('category')
        apart_mixed_leiden = adata_cell.obs['leiden'].iloc[np.where(adata_cell.obs['label'] != 'Mixed')].astype(
            'str').astype('category')
    else:
        apart_mixed_label = adata_cell.obs['label']
        apart_mixed_louvain = adata_cell.obs['louvain']
        apart_mixed_leiden = adata_cell.obs['leiden']

    score = pd.DataFrame({'ARI': [adjusted_rand_score(apart_mixed_louvain, apart_mixed_label),
                                  adjusted_rand_score(apart_mixed_leiden, apart_mixed_label)],
                          'NMI': [normalized_mutual_info_score(apart_mixed_louvain, apart_mixed_label),
                                  normalized_mutual_info_score(apart_mixed_leiden, apart_mixed_label)]},
                         index=['louvain', 'leiden'])
    print(score)

    fig3, ax = plt.subplots(figsize=(3, 3))
    sc.pl.umap(adata_cell, color='louvain', title='Louvain', s=20, ax=ax)
    if save_fig:
        plt.savefig(save_file + 'louvain.csv')

    fig4, ax = plt.subplots(figsize=(3, 3))
    sc.pl.umap(adata_cell, color='leiden', title='Leiden', s=20, ax=ax)
    if save_fig:
        plt.savefig(save_file + 'leiden.csv')

    # ========================================================================
    if plot_confusion:
        if data_name == 'GM12878vsHL':
            name = dict(zip(['GM12878', 'HL60', 'Mixed'], range(3)))
        elif data_name == 'GM12878vsHEK':
            name = dict(zip(['GM12878', 'HEK293T', 'Mixed'], range(3)))
        else:
            name = dict(zip(set(labels), range(cluster_num)))
        Cluster_louvain_mat = confusion_matrix(adata_cell.obs['louvain'].astype(int), [name[i] for i in labels])
        Cluster_leiden_mat = confusion_matrix(adata_cell.obs['leiden'].astype(int), [name[i] for i in labels])

        if data_name == 'GM12878vsHEK' or data_name == 'GM12878vsHL':
            Cluster_louvain_mat = pd.DataFrame(Cluster_louvain_mat[:-1, :],
                                               columns=list(name.keys()),
                                               index=range(cluster_num))
            Cluster_leiden_mat = pd.DataFrame(Cluster_leiden_mat[:-1, :],
                                              columns=list(name.keys()),
                                              index=range(cluster_num))
        else:
            Cluster_louvain_mat = pd.DataFrame(Cluster_louvain_mat, columns=list(name.keys()), index=range(cluster_num))
            Cluster_leiden_mat = pd.DataFrame(Cluster_leiden_mat, columns=list(name.keys()), index=range(cluster_num))

        Cluster_louvain_mat = confusion_sort(Cluster_louvain_mat)
        Cluster_leiden_mat = confusion_sort(Cluster_leiden_mat)

        sns.heatmap(Cluster_louvain_mat, annot=True, cmap="Blues", fmt='.20g', cbar=False)
        plt.xticks(rotation=30, horizontalalignment="center")
        plt.ylabel('Predicted')
        plt.xlabel('True')
        plt.title(data_name + '(louvain)', size=10)
        if save_fig:
            plt.savefig(save_file + 'confusion_louvain.csv')

        sns.heatmap(Cluster_leiden_mat, annot=True, cmap="Blues", fmt='.20g', cbar=False)
        plt.xticks(rotation=30, horizontalalignment="center")
        plt.ylabel('Predicted')
        plt.xlabel('True')
        plt.title(data_name + '(leiden)', size=10)
        if save_fig:
            plt.savefig(save_file + 'confusion_leiden.csv')
