import argparse
import os
import subprocess
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from utils import save


THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run original scGNN on h5ad data through the benchmark interface.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--data_path', type=str, required=True, help='Path to input h5ad file')
    parser.add_argument('--save_dir', type=str, default='./results', help='Directory to save results')
    parser.add_argument('--n_clusters', type=int, default=0, help='Cluster count for fixed-k clustering modes')
    parser.add_argument('--dataset_name', type=str, default=None, help='Optional dataset alias used for working files')

    parser.add_argument('--seed', type=int, default=1, help='Random seed passed to original scGNN')
    parser.add_argument('--gpu', type=int, default=0, help='GPU device ID')
    parser.add_argument('--no_cuda', action='store_true', help='Disable CUDA')
    parser.add_argument('--cores_usage', type=str, default='1', help='Number of CPU cores used by original scGNN')

    parser.add_argument('--gene_selectnum', type=int, default=2000, help='Number of genes selected in original preprocessing')
    parser.add_argument('--cell_ratio', type=float, default=0.99, help='Cell filter ratio in original preprocessing')
    parser.add_argument('--gene_ratio', type=float, default=0.99, help='Gene filter ratio in original preprocessing')
    parser.add_argument('--transform', type=str, default='log', help='Transform mode in original preprocessing')

    parser.add_argument('--model', type=str, default='AE', choices=['AE', 'VAE'], help='Original scGNN autoencoder type')
    parser.add_argument('--batch_size', type=int, default=12800, help='Original scGNN batch size')
    parser.add_argument('--Regu_epochs', type=int, default=500, help='Original scGNN pretrain epochs')
    parser.add_argument('--EM_epochs', type=int, default=200, help='Original scGNN EM epochs')
    parser.add_argument('--EM_iteration', type=int, default=10, help='Original scGNN EM iterations')
    parser.add_argument('--cluster_epochs', type=int, default=200, help='Original scGNN cluster AE epochs')
    parser.add_argument('--quickmode', action='store_true', help='Use original scGNN quickmode')

    parser.add_argument('--regulized_type', type=str, default='noregu', choices=['noregu', 'LTMG', 'LTMG01'], help='Original scGNN regularization type')
    parser.add_argument('--gammaPara', type=float, default=0.1, help='Original scGNN gammaPara')
    parser.add_argument('--alphaRegularizePara', type=float, default=0.9, help='Original scGNN alphaRegularizePara')
    parser.add_argument('--reduction', type=str, default='sum', choices=['sum', 'mean'], help='Original scGNN loss reduction')

    parser.add_argument('--k', type=int, default=10, help='KNN k in original scGNN graph building')
    parser.add_argument('--knn_distance', type=str, default='euclidean', help='Distance metric in original scGNN graph building')
    parser.add_argument('--prunetype', type=str, default='KNNgraphStatsSingleThread', help='Graph pruning type in original scGNN')

    parser.add_argument('--clustering_method', type=str, default='LouvainK', help='Original scGNN clustering method')
    parser.add_argument('--maxClusterNumber', type=int, default=30, help='Original scGNN maxClusterNumber')
    parser.add_argument('--minMemberinCluster', type=int, default=5, help='Original scGNN minMemberinCluster')
    parser.add_argument('--resolution', type=str, default='auto', help='Original scGNN Louvain resolution mode')

    parser.add_argument('--useGAEembedding', action='store_true', help='Use original scGNN GAE embedding')
    parser.add_argument('--useBothembedding', action='store_true', help='Use both feature and GAE embeddings')
    parser.add_argument('--GAEepochs', type=int, default=200, help='Original scGNN GAE epochs')
    parser.add_argument('--GAEhidden1', type=int, default=32, help='Original scGNN GAE hidden1 size')
    parser.add_argument('--GAEhidden2', type=int, default=16, help='Original scGNN GAE hidden2 size')
    parser.add_argument('--GAElr', type=float, default=0.01, help='Original scGNN GAE learning rate')
    parser.add_argument('--GAEdropout', type=float, default=0.0, help='Original scGNN GAE dropout')
    parser.add_argument('--GAEmodel', type=str, default='gcn_vae', help='Original scGNN GAE model')
    parser.add_argument('--GAElr_dw', type=float, default=0.001, help='Original scGNN GAE regularization learning rate')

    return parser.parse_args()


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def decode_vector(values):
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode('utf-8'))
        elif hasattr(value, 'decode'):
            decoded.append(value.decode('utf-8'))
        else:
            decoded.append(str(value))
    return decoded


def read_csr(group):
    shape = tuple(group.attrs['shape'])
    return sp.csr_matrix((group['data'][...], group['indices'][...], group['indptr'][...]), shape=shape)


def load_labels_from_h5ad(data_path, label_col='Celltype'):
    with h5py.File(data_path, 'r') as handle:
        obs = handle['obs']
        cell_ids = decode_vector(obs['_index'][...])
        if label_col not in obs:
            raise KeyError(f'Label column {label_col!r} not found in h5ad obs. Available columns: {list(obs.keys())}')
        labels = decode_vector(obs[label_col][...])
    return dict(zip(cell_ids, labels)), len(set(labels))


def export_h5ad_to_10x(data_path, output_root):
    with h5py.File(data_path, 'r') as handle:
        matrix = read_csr(handle['X']).transpose().tocoo().astype(np.int32)
        cell_ids = decode_vector(handle['obs']['_index'][...])
        gene_ids = decode_vector(handle['var']['_index'][...])

    ensure_dir(output_root)

    with open(os.path.join(output_root, 'barcodes.tsv'), 'w', encoding='utf-8') as file_obj:
        file_obj.write('\n'.join(cell_ids) + '\n')

    with open(os.path.join(output_root, 'features.tsv'), 'w', encoding='utf-8') as file_obj:
        file_obj.write('\n'.join(gene_ids) + '\n')

    order = np.lexsort((matrix.row, matrix.col))
    matrix_path = os.path.join(output_root, 'matrix.mtx')
    with open(matrix_path, 'w', encoding='utf-8') as file_obj:
        file_obj.write('%%MatrixMarket matrix coordinate integer general\n')
        file_obj.write(f'{matrix.shape[0]} {matrix.shape[1]} {matrix.nnz}\n')
        for idx in order:
            file_obj.write(f'{matrix.row[idx] + 1} {matrix.col[idx] + 1} {int(matrix.data[idx])}\n')


def run_command(command, cwd):
    print('Running:', ' '.join(command))
    subprocess.run(command, cwd=cwd, check=True)


def build_dataset_alias(args):
    if args.dataset_name:
        return args.dataset_name
    base = os.path.splitext(os.path.basename(args.data_path))[0]
    return f'{base}_scGNN'


def original_scgnn_paths(args, dataset_alias):
    save_dir = os.path.abspath(args.save_dir)
    work_root = ensure_dir(os.path.join(save_dir, 'workdir'))
    dataset_dir = ensure_dir(os.path.join(work_root, dataset_alias))
    original_output_dir = ensure_dir(os.path.join(save_dir, 'original_scgnn'))
    return work_root, dataset_dir, original_output_dir


def maybe_check_ltmg(args):
    if args.regulized_type == 'noregu':
        return
    try:
        import rpy2  # noqa: F401
    except Exception as exc:
        raise RuntimeError(
            'Original scGNN LTMG mode requires rpy2 and the R package scGNNLTMG, which are not available in this environment.'
        ) from exc


def preprocess_with_original_scgnn(args, dataset_alias, work_root):
    preprocess_command = [
        sys.executable,
        os.path.join(THIS_DIR, 'PreprocessingscGNN.py'),
        '--datasetName', dataset_alias,
        '--datasetDir', work_root + os.sep,
        '--LTMGDir', work_root + os.sep,
        '--filetype', '10X',
        '--geneSelectnum', str(args.gene_selectnum),
        '--cellRatio', str(args.cell_ratio),
        '--geneRatio', str(args.gene_ratio),
        '--transform', args.transform,
    ]
    if args.regulized_type != 'noregu':
        preprocess_command.append('--inferLTMGTag')
    run_command(preprocess_command, cwd=THIS_DIR)


def train_with_original_scgnn(args, dataset_alias, work_root, original_output_dir):
    command = [
        sys.executable,
        os.path.join(THIS_DIR, 'scGNN.py'),
        '--datasetName', dataset_alias,
        '--datasetDir', work_root + os.sep,
        '--outputDir', original_output_dir + os.sep,
        '--seed', str(args.seed),
        '--model', args.model,
        '--batch-size', str(args.batch_size),
        '--Regu-epochs', str(args.Regu_epochs),
        '--EM-epochs', str(args.EM_epochs),
        '--EM-iteration', str(args.EM_iteration),
        '--cluster-epochs', str(args.cluster_epochs),
        '--regulized-type', args.regulized_type,
        '--gammaPara', str(args.gammaPara),
        '--alphaRegularizePara', str(args.alphaRegularizePara),
        '--reduction', args.reduction,
        '--k', str(args.k),
        '--knn-distance', args.knn_distance,
        '--prunetype', args.prunetype,
        '--coresUsage', str(args.cores_usage),
        '--clustering-method', args.clustering_method,
        '--maxClusterNumber', str(args.maxClusterNumber),
        '--minMemberinCluster', str(args.minMemberinCluster),
        '--resolution', args.resolution,
        '--GAEepochs', str(args.GAEepochs),
        '--GAEhidden1', str(args.GAEhidden1),
        '--GAEhidden2', str(args.GAEhidden2),
        '--GAElr', str(args.GAElr),
        '--GAEdropout', str(args.GAEdropout),
        '--GAEmodel', args.GAEmodel,
        '--GAElr_dw', str(args.GAElr_dw),
    ]
    if args.n_clusters > 0:
        command.extend(['--n-clusters', str(args.n_clusters)])
    if args.no_cuda:
        command.append('--no-cuda')
    if args.quickmode:
        command.append('--quickmode')
    if args.useGAEembedding:
        command.append('--useGAEembedding')
    if args.useBothembedding:
        command.append('--useBothembedding')

    run_command(command, cwd=THIS_DIR)


def collect_and_save_results(data_path, dataset_alias, save_dir, original_output_dir):
    label_map, inferred_clusters = load_labels_from_h5ad(data_path)
    embedding_path = os.path.join(original_output_dir, f'{dataset_alias}_embedding.csv')
    result_path = os.path.join(original_output_dir, f'{dataset_alias}_results.txt')

    if not os.path.exists(embedding_path):
        raise FileNotFoundError(f'Expected original scGNN embedding file not found: {embedding_path}')
    if not os.path.exists(result_path):
        raise FileNotFoundError(f'Expected original scGNN result file not found: {result_path}')

    embedding_df = pd.read_csv(embedding_path, index_col=0)
    result_df = pd.read_csv(result_path, index_col=0)

    cell_ids = embedding_df.index.tolist()
    missing = [cell_id for cell_id in cell_ids if cell_id not in label_map]
    if missing:
        raise KeyError(f'{len(missing)} output cells are missing from h5ad labels, first example: {missing[0]}')

    y_true = np.array([label_map[cell_id] for cell_id in cell_ids])
    y_true = pd.factorize(y_true)[0]
    y_pred = pd.factorize(result_df.iloc[:, 0].astype(str))[0]
    embedding = embedding_df.to_numpy()

    save(save_dir, y_true, y_pred, 0, embedding)
    print(f'Original scGNN output clusters: {len(np.unique(y_pred))}')
    print(f'H5AD reference cell types: {inferred_clusters}')


def main():
    args = parse_args()
    maybe_check_ltmg(args)

    save_dir = ensure_dir(os.path.abspath(args.save_dir))
    dataset_alias = build_dataset_alias(args)
    work_root, dataset_dir, original_output_dir = original_scgnn_paths(args, dataset_alias)

    print('Exporting h5ad to 10X for original scGNN preprocessing...')
    export_h5ad_to_10x(args.data_path, dataset_dir)

    print('Running original scGNN preprocessing...')
    preprocess_with_original_scgnn(args, dataset_alias, work_root)

    print('Running original scGNN training and clustering...')
    train_with_original_scgnn(args, dataset_alias, work_root, original_output_dir)

    print('Collecting original scGNN outputs into benchmark format...')
    collect_and_save_results(args.data_path, dataset_alias, save_dir, original_output_dir)

    print(f'Finished. Benchmark outputs are in: {save_dir}')
    print(f'Original scGNN raw outputs are in: {original_output_dir}')


if __name__ == '__main__':
    main()
