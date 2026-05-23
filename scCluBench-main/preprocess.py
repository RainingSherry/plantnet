from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import pickle, os, numbers

import h5py
import numpy as np
import scipy as sp
import pandas as pd
import scanpy as sc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import scale
import utils as utils


# TODO: Fix this
class AnnSequence:
    def __init__(self, matrix, batch_size, sf=None):
        self.matrix = matrix
        if sf is None:
            self.size_factors = np.ones((self.matrix.shape[0], 1),
                                        dtype=np.float32)
        else:
            self.size_factors = sf
        self.batch_size = batch_size

    def __len__(self):
        return len(self.matrix) // self.batch_size

    def __getitem__(self, idx):
        batch = self.matrix[idx * self.batch_size:(idx + 1) * self.batch_size]
        batch_sf = self.size_factors[idx * self.batch_size:(idx + 1) * self.batch_size]

        # return an (X, Y) pair
        return {'count': batch, 'size_factors': batch_sf}, batch


def read_dataset(adata, transpose=False, test_split=False, copy=False):
    if isinstance(adata, sc.AnnData):
        if copy:
            adata = adata.copy()
    elif isinstance(adata, str):
        adata = sc.read(adata)
    else:
        raise NotImplementedError

    norm_error = 'Make sure that the dataset (adata.X) contains unnormalized count data.'
    assert 'n_count' not in adata.obs, norm_error

    if adata.X.size < 50e6:  # check if adata.X is integer only if array is small
        if sp.sparse.issparse(adata.X):
            assert (adata.X.astype(float) != adata.X).nnz == 0, norm_error
        else:
            assert np.all(adata.X.astype(float) == adata.X), norm_error

    if transpose: adata = adata.transpose()

    if test_split:
        train_idx, test_idx = train_test_split(np.arange(adata.n_obs), test_size=0.1, random_state=42)
        spl = pd.Series(['train'] * adata.n_obs)
        spl.iloc[test_idx] = 'test'
        adata.obs['DCA_split'] = spl.values
    else:
        adata.obs['DCA_split'] = 'train'

    adata.obs['DCA_split'] = adata.obs['DCA_split'].astype('category')
    print('### Autoencoder: Successfully preprocessed {} genes and {} cells.'.format(adata.n_vars, adata.n_obs))

    return adata


def normalize(adata, filter_min_counts=True, size_factors=True, normalize_input=True, logtrans_input=True):
    if filter_min_counts:
        sc.pp.filter_genes(adata, min_counts=1)
        sc.pp.filter_cells(adata, min_counts=1)

    if size_factors or normalize_input or logtrans_input:
        adata.raw = adata.copy()
    else:
        adata.raw = adata

    if size_factors:
        sc.pp.normalize_per_cell(adata)
        adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    else:
        adata.obs['size_factors'] = 1.0

    if logtrans_input:
        sc.pp.log1p(adata)

    if normalize_input:
        sc.pp.scale(adata)

    return adata


def read_genelist(filename):
    genelist = list(set(open(filename, 'rt').read().strip().split('\n')))
    assert len(genelist) > 0, 'No genes detected in genelist file'
    print('### Autoencoder: Subset of {} genes will be denoised.'.format(len(genelist)))

    return genelist


def write_text_matrix(matrix, filename, rownames=None, colnames=None, transpose=False):
    if transpose:
        matrix = matrix.T
        rownames, colnames = colnames, rownames

    pd.DataFrame(matrix, index=rownames, columns=colnames).to_csv(filename,
                                                                  sep='\t',
                                                                  index=(rownames is not None),
                                                                  header=(colnames is not None),
                                                                  float_format='%.6f')

def read_pickle(inputfile):
    return pickle.load(open(inputfile, "rb"))

# new
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score, fowlkes_mallows_score, v_measure_score, silhouette_score, accuracy_score
from sklearn.metrics.cluster import homogeneity_score, completeness_score
import numpy as np
import scanpy as sc
import scipy.sparse as sp

def check_normalization(adata):
    """
    检查 AnnData 是否已经归一化、log1p 变换和标准化
    Returns:
        - is_norm: 是否归一化
        - is_log1p: 是否 log1p 变换
        - is_scaled: 是否标准化
    """
    X = adata.X

    # **解决稀疏矩阵问题**
    if sp.issparse(X):  
        X = X.toarray()  # 转换为 NumPy 数组，避免 `ValueError`
    
    is_norm = not np.all(X.astype(int) == X)  # 如果全是整数，则未归一化
    is_log1p = X.max() < 20  # log1p 后最大值一般 <20
    is_scaled = np.isclose(X.mean(), 0, atol=0.1) and np.isclose(X.std(), 1, atol=0.1)  # 均值0，方差1

    return is_norm, is_log1p, is_scaled



def normalize_sc(adata, size_factors=True, filter_min_counts=True, logtrans_input=True, normalize_input=True):
    """
    归一化、log1p 变换和标准化 scRNA-seq 数据，并保存原始数据到 adata.raw
    """
    # Step 1: 确保 adata.raw 存储原始未归一化数据
    if adata.raw is None:
        adata.raw = adata.copy()

    # # Step 7: 过滤低表达基因
    # if filter_min_counts:
    #     sc.pp.filter_genes(adata, min_counts=3)
        
    if filter_min_counts:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes = 1000, subset=True)

    # Step 2: 检查当前数据状态
    is_norm, is_log1p, is_scaled = check_normalization(adata)
    print(f"是否归一化: {is_norm}, 是否 log1p 变换: {is_log1p}, 是否标准化: {is_scaled}")

    # Step 3: 如果未归一化，则进行归一化
    if not is_norm:
        print("数据未归一化，进行 normalize_per_cell 处理")
        sc.pp.normalize_per_cell(adata)  
        is_norm = True  # 更新状态

    # # Step 4: 计算 size factor
    # if size_factors:
    #     adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    # Step 4: 计算 size factor
    if size_factors:
        count_column = 'n_counts' if 'n_counts' in adata.obs.columns else 'total_counts'
        adata.obs['size_factors'] = adata.obs[count_column] / np.median(adata.obs[count_column])

    
    # Step 5: 如果未 log1p 变换，则进行 log1p 变换
    if logtrans_input and not is_log1p:
        print("数据未 log1p 变换，进行 log1p 处理")
        sc.pp.log1p(adata)
        is_log1p = True  # 更新状态

    # Step 6: 如果未标准化，则进行标准化
    if normalize_input and not is_scaled:
        print("数据未标准化，进行 scale 处理")
        sc.pp.scale(adata)
        is_scaled = True  # 更新状态

    # print(adata)
    # print(adata.var.highly_variable.index)
    # print(adata.var.highly_variable)

    return adata



def prepare_data_for_model(file_path, size_factors=True, filter_min_counts=True, logtrans_input=True, normalize_input=True):
    """
    This function reads the data from a file, normalizes, and prepares the data for modeling.

    Parameters:
    - file_path: Path to the input file.
    - size_factors: Whether to use size factors for normalization.
    - highly_genes: Whether to use only highly variable genes.
    - logtrans_input: Whether to perform log-transformation.
    - normalize_input: Whether to scale the data.

    Returns:
    - X: DataFrame of normalized gene expression data.
    - Y: Series of cell types.
    - sf: Size factors.
    """

    # Read data
    data = sc.read_h5ad(file_path)

    # Normalize the data
    data = normalize_sc(data, size_factors=size_factors, filter_min_counts=filter_min_counts, 
                         logtrans_input=logtrans_input, normalize_input=normalize_input)

    # Prepare the data for the model
    X = data.to_df()  # Get data as DataFrame
    Y = data.obs['cell_type']  # Get cell type labels
    sf = data.obs['size_factors']  # Get size factors
    

    return X, Y, sf, data

def read_clean(data):
    assert isinstance(data, np.ndarray)
    if data.dtype.type is np.bytes_:
        data = utils.decode(data)
    if data.size == 1:
        data = data.flat[0]
    return data


def dict_from_group(group):
    assert isinstance(group, h5py.Group)
    d = utils.dotdict()
    for key in group:
        if isinstance(group[key], h5py.Group):
            value = dict_from_group(group[key])
        else:
            value = read_clean(group[key][...])
        d[key] = value
    return d


def read_data(filename, sparsify = False, skip_exprs = False):
    with h5py.File(filename, "r") as f:
        obs = pd.DataFrame(dict_from_group(f["obs"]), index = utils.decode(f["obs_names"][...]))
        var = pd.DataFrame(dict_from_group(f["var"]), index = utils.decode(f["var_names"][...]))
        uns = dict_from_group(f["uns"])
        if not skip_exprs:
            exprs_handle = f["exprs"]
            if isinstance(exprs_handle, h5py.Group):
                mat = sp.csr_matrix((exprs_handle["data"][...], exprs_handle["indices"][...],
                                               exprs_handle["indptr"][...]), shape = exprs_handle["shape"][...])
            else:
                mat = exprs_handle[...].astype(np.float32)
                if sparsify:
                    mat = sp.csr_matrix(mat)
        else:
            mat = sp.csr_matrix((obs.shape[0], var.shape[0]))
    return mat, obs, var, uns


def prepro(data_path):
    data_path = data_path

    mat, obs, var, uns = read_data(data_path, sparsify=False, skip_exprs=False)
    if isinstance(mat, np.ndarray):
        X = np.array(mat)
    else:
        X = np.array(mat.toarray())
    cell_name = np.array(obs["cell_type1"])
    cell_type, cell_label = np.unique(cell_name, return_inverse=True)
    return X, cell_label, var, cell_name

def normalize(adata, copy=True, highly_genes = None, filter_min_counts=True, size_factors=True, normalize_input=True, logtrans_input=True):
    if isinstance(adata, sc.AnnData):
        if copy:
            adata = adata.copy()
    elif isinstance(adata, str):
        adata = sc.read(adata)
    else:
        raise NotImplementedError
    norm_error = 'Make sure that the dataset (adata.X) contains unnormalized count data.'
    assert 'n_count' not in adata.obs, norm_error
    if adata.X.size < 50e6: # check if adata.X is integer only if array is small
        if sp.sparse.issparse(adata.X):
            assert (adata.X.astype(int) != adata.X).nnz == 0, norm_error
        else:
            assert np.all(adata.X.astype(int) == adata.X), norm_error

    if filter_min_counts:
        sc.pp.filter_genes(adata, min_counts=1)
        sc.pp.filter_cells(adata, min_counts=1)
    if size_factors or normalize_input or logtrans_input:
        adata.raw = adata.copy()
    else:
        adata.raw = adata
    if size_factors:
        sc.pp.normalize_per_cell(adata)
        adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    else:
        adata.obs['size_factors'] = 1.0
    if logtrans_input:
        sc.pp.log1p(adata)
    if highly_genes != None:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes = highly_genes, subset=True)
    if normalize_input:
        sc.pp.scale(adata)
    return adata