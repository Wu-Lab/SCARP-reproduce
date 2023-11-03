import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

merge_thre_louvain_NMI_score = pd.read_table('./Results/merge_thre sensitivity louvain_NMI.txt', sep='\t', index_col=0)
m_louvain_NMI_score = pd.read_table('./Results/m sensitivity louvain_NMI.txt', sep='\t', index_col=0)
beta_louvain_NMI_score = pd.read_table('./Results/beta sensitivity louvain_NMI.txt', sep='\t', index_col=0)

merge_thre_louvain_ARI_score = pd.read_table('./Results/merge_thre sensitivity louvain_ARI.txt', sep='\t', index_col=0)
m_louvain_ARI_score = pd.read_table('./Results/m sensitivity louvain_ARI.txt', sep='\t', index_col=0)
beta_louvain_ARI_score = pd.read_table('./Results/beta sensitivity louvain_ARI.txt', sep='\t', index_col=0)

louvain_NMI_score = pd.DataFrame({
    'NMI scores': list(np.array(merge_thre_louvain_NMI_score).astype('float').flatten()) +
                  list(np.array(beta_louvain_NMI_score).astype('float').flatten()) +
                  list(np.array(m_louvain_NMI_score).astype('float').flatten()),
    'datasets': list(np.repeat(list(merge_thre_louvain_NMI_score.index), merge_thre_louvain_NMI_score.shape[1])) +
                list(np.repeat(list(beta_louvain_NMI_score.index), beta_louvain_NMI_score.shape[1])) +
                list(np.repeat(list(m_louvain_NMI_score.index), m_louvain_NMI_score.shape[1])),
    'parameters': list(np.repeat('γ', 9 * merge_thre_louvain_NMI_score.shape[1])) +
                  list(np.repeat('β', 9 * beta_louvain_NMI_score.shape[1])) +
                  list(np.repeat('m', 9 * m_louvain_NMI_score.shape[1]))})
louvain_NMI_score['parameters'] = louvain_NMI_score['parameters'].astype('category')

louvain_ARI_score = pd.DataFrame({
    'ARI scores': list(np.array(merge_thre_louvain_ARI_score).astype('float').flatten()) +
                  list(np.array(beta_louvain_ARI_score).astype('float').flatten()) +
                  list(np.array(m_louvain_ARI_score).astype('float').flatten()),
    'datasets': list(np.repeat(list(merge_thre_louvain_ARI_score.index), merge_thre_louvain_ARI_score.shape[1])) +
                list(np.repeat(list(beta_louvain_ARI_score.index), beta_louvain_ARI_score.shape[1])) +
                list(np.repeat(list(m_louvain_ARI_score.index), m_louvain_ARI_score.shape[1])),
    'parameters': list(np.repeat('γ', 9 * merge_thre_louvain_ARI_score.shape[1])) +
                  list(np.repeat('β', 9 * beta_louvain_ARI_score.shape[1])) +
                  list(np.repeat('m', 9 * m_louvain_ARI_score.shape[1]))})
louvain_ARI_score['parameters'] = louvain_ARI_score['parameters'].astype('category')

plot_data = louvain_ARI_score.copy()
plot_data['NMI scores'] = louvain_NMI_score['NMI scores']
plot_data['Average of NMI and ARI scores'] = (plot_data['ARI scores']+plot_data['NMI scores'])/2

sns.set_style('whitegrid')
plt.figure(figsize=(7, 2.7))
sns.barplot(x='datasets', y='Average of NMI and ARI scores',
            hue='parameters', data=plot_data,
            palette="Accent")
plt.title('Sensitivity to parameter selection')
plt.legend(bbox_to_anchor=(1.01, 0.3), loc=3, borderaxespad=0)
plt.xticks(rotation=35)
plt.subplots_adjust(bottom=0.4, left=0.1, right=0.8)
plt.show()
plt.savefig('./Results/para_mean.svg')

sns.set_style('whitegrid')
plt.figure(figsize=(7, 2.7))
sns.barplot(x='datasets', y='NMI scores',
            hue='parameters', data=plot_data,
            palette="Accent")
plt.title('Sensitivity to parameter selection')
plt.legend(bbox_to_anchor=(1.01, 0.3), loc=3, borderaxespad=0)
plt.xticks(rotation=35)
plt.subplots_adjust(bottom=0.4, left=0.1, right=0.8)
plt.show()
plt.savefig('./Results/para_NMI.svg')


sns.set_style('whitegrid')
plt.figure(figsize=(7, 2.7))
sns.barplot(x='datasets', y='ARI scores',
            hue='parameters', data=plot_data,
            palette="Accent")
plt.title('Sensitivity to parameter selection')
plt.legend(bbox_to_anchor=(1.01, 0.3), loc=3, borderaxespad=0)
plt.xticks(rotation=35)
plt.subplots_adjust(bottom=0.4, left=0.1, right=0.8)
plt.show()
plt.savefig('./Results/para_ARI.svg')
