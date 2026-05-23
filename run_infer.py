#!/usr/bin/env python3
import sys, os, subprocess
os.environ["TMPDIR"] = "/data/tmp"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

result = subprocess.run(
    ["/data/luolie/conda/envs/scclubench-main/bin/python",
     "/home/luolie/biopipeline/dimension-reduction/plantnet/infer_graphdiffusion_ckpt.py",
     "--data-path", "/home/luolie/biopipeline/dimension-reduction/plantnet/data/SRP182008.h5ad",
     "--result-dir", "/home/luolie/biopipeline/dimension-reduction/plantnet/results/cursor_Doloris_GraphDiffusion/SRP182008",
     "--config-path", "/home/luolie/biopipeline/dimension-reduction/plantnet/results/cursor_Doloris_GraphDiffusion/SRP182008/config.json",
     "--checkpoint-path", "/home/luolie/biopipeline/dimension-reduction/plantnet/results/cursor_Doloris_GraphDiffusion/SRP182008/best_model.pt",
     "--n-clusters", "15",
     "--gpu", "1"],
    capture_output=True, text=True
)

with open("/data/tmp/gd_out.txt", "w") as f:
    f.write("STDOUT:\n" + result.stdout + "\nSTDERR:\n" + result.stderr + "\nRC=" + str(result.returncode) + "\n")
print("Done, wrote to /data/tmp/gd_out.txt")
