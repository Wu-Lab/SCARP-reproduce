import sys
sys.path.append('../Scarp/')

from Scarp.downstream import *
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")
np.random.seed(1)
random_state = 1
beta = 5000
m = 1.5
p = 1


# ===========================================================================================
#                                        import data
# ===========================================================================================
Data_name = ['Forebrain', 'GM12878vsHL', 'Leukemia',
             'Breast_Tumor', 'GM12878vsHEK', 'InSilico',
             'Sox10KD_filter_var0.6', 'blood2K_filter60', 'Splenocyte']

merge_thre_value = np.append(np.array([100, 500]), np.arange(1000, 5001, 1000))

louvain_NMI_score = pd.DataFrame(index=Data_name, columns=merge_thre_value)
louvain_ARI_score = pd.DataFrame(index=Data_name, columns=merge_thre_value)
Kept_component = []

for d in range(len(Data_name)):
    data_name = Data_name[d]
    # k = Kept_comp[d]
    print('\n++++++++++++++++++++ %s ++++++++++++++++++++++++' % data_name)

    # ===========================================================================================
    data = sc.read_h5ad('../Exp1_Benchmark/Processed Data/' + data_name + '.h5ad')
    Cells = data.obs.index
    Cells_num, Peaks_num = data.X.shape
    N = Cells_num + Peaks_num
    labels = data.obs['celltype'].astype('category')
    if data_name == 'GM12878vsHEK' or data_name == 'GM12878vsHL':
        cluster_num = 2
    else:
        cluster_num = np.unique(labels).shape[0]
    print('Number of Peaks:', Peaks_num)
    print('Number of Cells:', Cells_num)
    print('Number of labels: ', cluster_num)

    for j in range(merge_thre_value.shape[0]):
        merge_thre = merge_thre_value[j]
        t, diffusion_mat = SCARP(data=data,
                                 m=m,
                                 merge_thre=merge_thre,
                                 beta=beta
                                 )

        if data_name == 'Breast_Tumor':
            k = 5
        elif data_name == 'GM12878vsHL':
            k = 5
        elif data_name == 'GM12878vsHEK':
            k = 5
        else:
            if Peaks_num > 50000:
                k = std_plot(data=diffusion_mat,
                             title=data_name,
                             max_k=100,
                             plot_std=False)
            else:
                k = std_plot(data=diffusion_mat,
                             title=data_name,
                             max_k=50,
                             plot_std=False)
        Kept_component.append(k)

        SCARP_score = compute_score(data_name=data_name,
                                    diffusion_mat=diffusion_mat,
                                    cluster_num=cluster_num,
                                    labels=labels,
                                    Cells=Cells,
                                    kept_comp=k,
                                    random_state=1
                                    )
        louvain_NMI_score[merge_thre][data_name] = SCARP_score['NMI']['louvain']
        louvain_ARI_score[merge_thre][data_name] = SCARP_score['ARI']['louvain']



# =========================================================================================
#                                   PLOT
# =========================================================================================
louvain_NMI_score = louvain_NMI_score.astype('float')
louvain_NMI_score1 = pd.DataFrame(scipy.stats.zscore(louvain_NMI_score, axis=1),
                                  index=Data_name, columns=merge_thre_value)
louvain_NMI_score1 = louvain_NMI_score1.fillna(0)

louvain_ARI_score = louvain_ARI_score.astype('float')
louvain_ARI_score1 = pd.DataFrame(scipy.stats.zscore(louvain_ARI_score, axis=1),
                                  index=Data_name, columns=merge_thre_value)
louvain_ARI_score1 = louvain_ARI_score1.fillna(0)

# ==============================
sns.set_style("darkgrid")
plt.figure(1, figsize=(10, 4))
sns.heatmap(louvain_NMI_score1, annot=louvain_NMI_score,
            linewidths=0.05, fmt=".3f",
            cbar=True,
            cmap='Blues',
            center=0,
            annot_kws={'color': 'black'})
plt.xticks(rotation=30, horizontalalignment="center")
plt.title('(Clustering: Louvain; Evaluation: NMI)')
plt.xlabel('merge threshold')
plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.17)
plt.show()
plt.savefig('./Results/merge_thre sensitivity louvain_NMI.svg', bbox_inches='tight')
louvain_NMI_score.to_csv('./Results/merge_thre sensitivity louvain_NMI.txt', sep='\t')

plt.figure(2, figsize=(10, 4))
sns.heatmap(louvain_ARI_score1, annot=louvain_ARI_score,
            linewidths=0.05, fmt=".3f",
            cbar=True,
            cmap='Blues',
            center=0,
            annot_kws={'color': 'black'})
plt.xticks(rotation=30, horizontalalignment="center")
plt.title('(Clustering: Louvain; Evaluation: ARI)')
plt.xlabel('merge threshold')
plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.17)
plt.show()
plt.savefig('./Results/merge_thre sensitivity louvain_ARI.svg', bbox_inches='tight')
louvain_ARI_score.to_csv('./Results/merge_thre sensitivity louvain_ARI.txt', sep='\t')

np.savetxt('./Results/merge_thre sensitivity Kept_component.txt', Kept_component)