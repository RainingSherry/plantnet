from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str]):
    print('running:', ' '.join(cmd))
    # Inherit PYTHONPATH so child scripts can find sibling modules
    env = os.environ.copy()
    child_path = str(Path(__file__).resolve().parent.parent)
    if 'PYTHONPATH' in env:
        env['PYTHONPATH'] = f"{child_path}:{env['PYTHONPATH']}"
    else:
        env['PYTHONPATH'] = child_path
    subprocess.run(cmd, check=True, env=env)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-path', type=str, required=True)
    parser.add_argument('--output-dir', type=str, required=True)
    parser.add_argument('--n-hvg', type=int, default=2000)
    parser.add_argument('--latent-dim', type=int, default=64)
    parser.add_argument('--hidden-dim', type=int, default=256)
    parser.add_argument('--batch-size', type=int, default=256)
    parser.add_argument('--source-epochs', type=int, default=20)
    parser.add_argument('--target-epochs', type=int, default=20)
    parser.add_argument('--bridge-epochs', type=int, default=30)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--diffusion-steps', type=int, default=50)
    parser.add_argument('--recon-weight', type=float, default=1.0)
    parser.add_argument('--prior-weight', type=float, default=1e-3)
    parser.add_argument('--zero-weight', type=float, default=0.25)
    parser.add_argument('--teacher-weight', type=float, default=1.0)
    parser.add_argument('--cluster-weight', type=float, default=1.0)
    parser.add_argument('--gaussian-weight', type=float, default=1e-3)
    parser.add_argument('--support-weight', type=float, default=0.5)
    parser.add_argument('--entropy-weight', type=float, default=1e-3)
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    common_source = [
        '--data-path', args.data_path,
        '--output-dir', str(output_dir),
        '--n-hvg', str(args.n_hvg),
        '--latent-dim', str(args.latent_dim),
        '--hidden-dim', str(args.hidden_dim),
        '--batch-size', str(args.batch_size),
        '--lr', str(args.lr),
        '--diffusion-steps', str(args.diffusion_steps),
        '--recon-weight', str(args.recon_weight),
        '--prior-weight', str(args.prior_weight),
        '--zero-weight', str(args.zero_weight),
    ]

    common_target = [
        '--data-path', args.data_path,
        '--output-dir', str(output_dir),
        '--n-hvg', str(args.n_hvg),
        '--latent-dim', str(args.latent_dim),
        '--hidden-dim', str(args.hidden_dim),
        '--batch-size', str(args.batch_size),
        '--lr', str(args.lr),
        '--diffusion-steps', str(args.diffusion_steps),
        '--recon-weight', str(args.recon_weight),
        '--prior-weight', str(args.prior_weight),
    ]

    common_bridge = [
        '--data-path', args.data_path,
        '--output-dir', str(output_dir),
        '--n-hvg', str(args.n_hvg),
        '--latent-dim', str(args.latent_dim),
        '--hidden-dim', str(args.hidden_dim),
        '--batch-size', str(args.batch_size),
        '--lr', str(args.lr),
        '--diffusion-steps', str(args.diffusion_steps),
        '--zero-weight', str(args.zero_weight),
        '--teacher-weight', str(args.teacher_weight),
        '--cluster-weight', str(args.cluster_weight),
        '--gaussian-weight', str(args.gaussian_weight),
        '--support-weight', str(args.support_weight),
        '--entropy-weight', str(args.entropy_weight),
    ]

    run([sys.executable, str(root / 'train_source.py'), *common_source, '--epochs', str(args.source_epochs)])
    run([sys.executable, str(root / 'train_target.py'), *common_target, '--epochs', str(args.target_epochs)])
    run([sys.executable, str(root / 'train_bridge_cluster.py'), *common_bridge, '--epochs', str(args.bridge_epochs)])


if __name__ == '__main__':
    main()
