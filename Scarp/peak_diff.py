import sys
sys.path.append('./Scarp/')
from downstream import *

def single_chromosome_diffusion_for_peaks_chain_analysis(groups, i, peak_to, sparse_matrix, Cells_num, beta, m):
    groups_i = groups.get_group(i)  # Current chromosome peak information
    peak_num_i = groups_i.shape[0]
    if '_' not in i:
        print('========== Chr' + str(i) + ': ' +
              str(peak_num_i) + ' Peaks ============')

    peak_from = peak_to
    peak_to = peak_from + peak_num_i
    #  =======Bipartite graph of the current chromosome===========
    sparse_i = sparse_matrix[:, peak_from:peak_to]
    sparse_i = np.array(sparse_i.todense())
    sparse_mat_i = np.vstack([np.hstack([np.zeros((Cells_num, Cells_num)), sparse_i]),
                              np.hstack([sparse_i.T, np.zeros((peak_num_i, peak_num_i))])])

    #  =======The isolated points of bipartite graph are removed first, and then added back later===========
    kept_index = np.where(sparse_mat_i.sum(1) == 0)[0]
    if kept_index.shape[0] != 0:
        sparse_mat_i = np.delete(sparse_mat_i, kept_index, axis=0)
        sparse_mat_i = np.delete(sparse_mat_i, kept_index, axis=1)
        # ==============================================================
        # # Padding peaks information
        peak_max = sparse_mat_i.max()
        Padding_mat = sparse_mat_i.copy()
        padding = create_peak_loc_subdiag_mat(
            groups_i, peak_max, beta)  # peaks*peaks matrix
        Padding_mat[(Cells_num - kept_index.shape[0]):, (Cells_num - kept_index.shape[0]):] = padding
        # ==============================================================
        diffusion_ = NR_F(Padding_mat, m)
        diffusion_ = np.array(diffusion_.todense())
        for kk in range(kept_index.shape[0]):
            diffusion_ = np.insert(diffusion_, kept_index[kk],
                                   np.zeros(diffusion_.shape[0]), axis=1)
            diffusion_ = np.insert(diffusion_, kept_index[kk],
                                   np.zeros(diffusion_.shape[1]), axis=0)
    else:
        # ==============================================================
        # Padding peaks information
        peak_max = sparse_mat_i.max()
        Padding_mat = sparse_mat_i.copy()
        padding = create_peak_loc_subdiag_mat(
            groups_i, peak_max, beta)  # peaks*peaks matrix
        Padding_mat[Cells_num:, Cells_num:] = padding
        # ==============================================================
        diffusion_ = NR_F(Padding_mat, m)
        diffusion_ = np.array(diffusion_.todense())

    return padding, diffusion_


# 从数据中提取某条染色体看连锁情况
def chain_analysis(data, chr_id, m, beta):
    Cells_num, Peaks_num = data.X.shape
    sparse_matrix = data.X  # sparse matrix
    sparse_matrix = (sparse_matrix > 0) * 1  # binary
    sparse_matrix = Mat_normalization(sparse_matrix)  # normalization

    # ==================================================================
    # Peaks information
    Peaks = data.var.index.tolist()
    Peaks_array = np.array([re.split('_|\W+', i) for i in Peaks])
    chri2i = dict(zip(['chr' + str(i) for i in range(1, 23)], range(1, 23)))
    chri2i['chrX'] = 23
    chri2i['chrY'] = 24
    Peaks_array[:, 0] = [chri2i[i] for i in Peaks_array[:, 0]]
    Peaks_df = pd.DataFrame(Peaks_array.astype(
        'float'), columns=['chr', 'from', 'to'])
    Peaks_df = Peaks_df.sort_values(by=['chr', 'from'])

    Peaks_df['chr'] = Peaks_df['chr'].astype('int').astype('str')
    chrs = Peaks_df['chr'].unique()
    groups = Peaks_df.groupby('chr')

    peak_to = [0]
    for i in chrs:
        peak_to.append(peak_to[-1] + groups.get_group(i).shape[0])
    peak_to = peak_to[0:-1]

    before_diffusion_PP, after_diffusion = single_chromosome_diffusion_for_peaks_chain_analysis(groups, str(chr_id),
                                                                                                peak_to[chr_id],
                                                                                                sparse_matrix,
                                                                                                Cells_num, beta, m)

    return before_diffusion_PP, after_diffusion, groups.get_group(str(chr_id))


def filter_peaks(adata, keep_peaks):
    adata_df = pd.DataFrame(adata.X.todense(),
                            index=adata.obs.index,
                            columns=adata.var.index)
    adata_new = sc.AnnData(adata_df[keep_peaks])
    adata_new.obs = adata.obs
    adata_new.var = adata.var.loc[keep_peaks]
    adata_new.X = csr_matrix(adata_new.X)
    return adata_new