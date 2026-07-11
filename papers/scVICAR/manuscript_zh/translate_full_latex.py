from __future__ import annotations

import json
import re
import shutil
import urllib.request
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[3]
EN = PROJECT / "papers/scVICAR/manuscript"
ZH = PROJECT / "papers/scVICAR/manuscript_zh/full_latex"
MODEL = "qwen3.6-27b-q4:latest"

FILES = [
    *(f"sections/{name}" for name in [
        "00_abstract.tex", "01_introduction.tex", "02_related_work.tex",
        "03_method.tex", "04_theory.tex", "05_experiments.tex",
        "06_results.tex", "07_downstream.tex", "08_discussion.tex",
        "08_availability.tex", "09_conclusion.tex",
    ]),
    *(f"generated/{name}" for name in [
        "baseline_results.tex", "confirmatory_results.tex", "downstream_results.tex",
        "full_label_sensitivity_results.tex", "leiden_results.tex", "stress_results.tex",
    ]),
]

TERM_RULES = """
固定术语：single-cell RNA sequencing=单细胞RNA测序；masked autoencoder=掩码自编码器；
graph-vicinal anchor recovery=图邻域锚点恢复；topology-informed affinity=拓扑信息亲和度；
cell-specific mixing coefficient=逐细胞混合系数；development benchmark=开发阶段基准；
frozen external comparison=冻结外部比较；matched-backbone=骨干匹配；
marker coherence=标记基因一致性；low-label transductive annotation=低标签量传导式注释。
方法名、数据集名、基因名、软件名、指标缩写、数学符号保持英文。
"""


def signature(text: str) -> dict[str, list[str]]:
    patterns = {
        "labels": r"\\label\{[^}]+\}",
        "refs": r"\\(?:ref|eqref|cref|Cref)\{[^}]+\}",
        "cites": r"\\cite\w*\{[^}]+\}",
        "inputs": r"\\input\{[^}]+\}",
        "graphics": r"\\includegraphics(?:\[[^]]*\])?\{[^}]+\}",
        "environments": r"\\(?:begin|end)\{[^}]+\}",
        "numbers": r"(?<![A-Za-z])[-+]?\d+(?:[,.]\d+)*(?:--\d+(?:\.\d+)?)?\\?%?",
    }
    return {key: sorted(re.findall(pattern, text)) for key, pattern in patterns.items()}


def validate(source: str, translated: str) -> list[str]:
    a, b = signature(source), signature(translated)
    return [key for key in a if a[key] != b[key]]


def generate(prompt: str) -> str:
    payload = json.dumps({
        "model": MODEL,
        "prompt": prompt,
        "stream": False,
        "think": False,
        "options": {"temperature": 0.05, "num_ctx": 32768},
    }).encode()
    request = urllib.request.Request(
        "http://127.0.0.1:11434/api/generate", data=payload,
        headers={"Content-Type": "application/json"}, method="POST",
    )
    with urllib.request.urlopen(request, timeout=900) as response:
        result = json.load(response)["response"].strip()
    result = re.sub(r"^```(?:latex|tex)?\s*", "", result)
    result = re.sub(r"\s*```$", "", result)
    return result.strip() + "\n"


def translate_file(relative: str) -> None:
    source_path = EN / relative
    target_path = ZH / relative
    source = source_path.read_text(encoding="utf-8")
    prompt = f"""你是计算生物学论文翻译编辑。把下面的英文 LaTeX 文件完整翻译成正式、直接、自然的学术中文。
严格要求：
1. 逐句翻译全部叙述文字，不概括、不删减、不增加论断。
2. 原样保留全部 LaTeX 命令、环境、公式、label、ref、cite、input、图片路径和注释。
3. 原样保留每个数字、正负号、百分比、区间、方法名和数据集名。
4. 翻译 section/subsection 标题、caption、表头和正文；数学模式内不翻译。
5. 不使用“不是X而是Y”模板，不添加解释或 Markdown 代码围栏。
{TERM_RULES}

源文件：
{source}"""
    translated = generate(prompt)
    problems = validate(source, translated)
    if problems:
        repair = f"""修复下面的中文 LaTeX，使它与英文源文件的结构签名完全一致。
不改变中文翻译措辞，只恢复遗漏或被改动的 LaTeX 项和数字。
不输出代码围栏。异常类别：{', '.join(problems)}

英文源：
{source}

待修复中文：
{translated}"""
        translated = generate(repair)
        problems = validate(source, translated)
    if problems:
        raise RuntimeError(f"{relative}: parity failure {problems}")
    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.write_text(translated, encoding="utf-8")
    print(f"PASS {relative}", flush=True)


def main() -> None:
    ZH.mkdir(parents=True, exist_ok=True)
    for relative in FILES:
        translate_file(relative)
    shutil.copy2(EN / "references.bib", ZH / "references.bib")
    shutil.copy2(EN / "results_macros.tex", ZH / "results_macros.tex")


if __name__ == "__main__":
    main()
