from __future__ import annotations

import re
import os
import subprocess
from pathlib import Path

import fitz


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[2]
SOURCE = ROOT / "scVICAR_中文阅读稿.md"
TEX = ROOT / "scVICAR_中文阅读稿.tex"
BUILD = ROOT / "build_latex"
PDF = ROOT / "scVICAR_中文阅读稿.pdf"
ASCII_PDF = ROOT / "scVICAR_zh_reading.pdf"
VECTOR_PDF = ROOT / "scVICAR_中文阅读稿_vector.pdf"
TECTONIC = PROJECT / ".codex_tex/bin/tectonic"


FIGURES = {
    "3. 图邻域锚点恢复方法": (
        "fig1_method.pdf",
        "scVICAR 方法框架。scVICAR-F 使用固定混合系数构造有界邻域视图；"
        "scVICAR-T 增加拓扑信息亲和度和逐细胞混合系数。两种变体均从邻域破坏视图恢复原始锚点。",
    ),
    "6.3 骨干匹配归因": (
        "fig2_confirmatory.pdf",
        "六个冻结数据集上的骨干匹配比较。图中给出相对 NoMix 的数据集层面 ARI 变化、"
        "变体分布、预设配对效应及共同训练预算下的运行时间。",
    ),
    "6.4 组件证据与负对照": (
        "fig3_components.pdf",
        "组件与负对照实验。RandomMix、固定邻域破坏、拓扑边权、逐细胞混合系数及完整模型"
        "使用相同骨干与训练预算。",
    ),
    "6.6 拓扑评分能够检测错误邻边风险": (
        "fig4_graph_stress.pdf",
        "三个数据集上的受控图污染实验。结果包括 ARI 退化曲线、scVICAR-T 与 scVICAR-F 的"
        "相对响应、边亲和度 AUROC 及逐细胞混合系数与邻域纯度的关联。",
    ),
    "7.4 结果": (
        "fig5_downstream.pdf",
        "六个冻结数据集上的下游效用。图中分别报告 marker recovery、marker-overlap annotation"
        " 和低标签量冻结表征线性探针。",
    ),
    "7.5 Human Pancreas 3 案例": (
        "fig6_pancreas_case.pdf",
        "Human Pancreas 3 预设案例。图中比较 NoMix、scVICAR-F 和 scVICAR-T 的 marker overlap、"
        "文献核验 marker panel 及类型到聚类和注释的对应关系。",
    ),
}


PREAMBLE = r"""\documentclass[journal]{IEEEtran}
\usepackage{amsmath,amssymb,amsthm,bm}
\usepackage{booktabs,multirow}
\usepackage{graphicx}
\usepackage{xcolor}
\usepackage{url}
\usepackage[hidelinks]{hyperref}
\usepackage{microtype}
\usepackage{placeins}
\usepackage[fontset=fandol]{ctex}
\graphicspath{{../figures/final/}}
\makeatletter
\setlength{\@fptop}{0pt}
\setlength{\@fpsep}{12pt}
\setlength{\@fpbot}{0pt plus 1fil}
\setlength{\@dblfptop}{0pt}
\setlength{\@dblfpsep}{12pt}
\setlength{\@dblfpbot}{0pt plus 1fil}
\makeatother
\renewcommand{\topfraction}{0.95}
\renewcommand{\floatpagefraction}{0.80}
\title{scVICAR：用于掩码单细胞表征学习的有界图邻域锚点恢复}
\author{匿名作者\thanks{用于双盲评审的中文阅读稿。投稿与引用以英文正式稿为准。}}
\begin{document}
\maketitle
"""


def escape_text(text: str) -> str:
    stash: dict[str, str] = {}

    def keep(pattern: str, prefix: str, repl=None) -> None:
        nonlocal text
        def sub(match: re.Match[str]) -> str:
            key = f"ZZ{prefix}{len(stash)}ZZ"
            stash[key] = repl(match) if repl else match.group(0)
            return key
        text = re.sub(pattern, sub, text)

    keep(r"\$[^$\n]+\$", "M")
    keep(r"`([^`]+)`", "C", lambda m: r"\texttt{" + latex_escape(m.group(1)) + "}")
    keep(r"\*\*([^*]+)\*\*", "B", lambda m: r"\textbf{" + latex_escape(m.group(1)) + "}")
    text = latex_escape(text)
    for key, value in stash.items():
        text = text.replace(key, value)
    return text


def latex_escape(text: str) -> str:
    return (text.replace("\\", r"\textbackslash{}")
            .replace("&", r"\&").replace("%", r"\%").replace("#", r"\#")
            .replace("_", r"\_").replace("{", r"\{").replace("}", r"\}"))


def figure_tex(name: str, caption: str, index: int) -> str:
    return rf"""
\begin{{figure*}}[t]
  \centering
  \includegraphics[width=\textwidth]{{{name}}}
  \caption{{{caption}}}
  \label{{fig:zh-{index}}}
\end{{figure*}}
"""


def convert(source: str) -> tuple[str, str, str]:
    source = re.sub(r"^# .*?\n", "", source, count=1)
    source = re.sub(r"^> .*?\n", "", source, count=1, flags=re.M)
    source = source.split("\n## 作者修改提示", 1)[0].strip()
    abstract_text, body = source.split("## 1. 引言", 1)
    abstract_text = abstract_text.replace("## 摘要", "", 1).strip()
    keyword_match = re.search(r"\*\*关键词：\*\*(.*)$", abstract_text, flags=re.M)
    keywords = keyword_match.group(1).strip() if keyword_match else ""
    if keyword_match:
        abstract_text = abstract_text[:keyword_match.start()].strip()
    abstract = "\n\n".join(escape_text(p.strip()) for p in abstract_text.split("\n\n") if p.strip())
    body = "## 1. 引言\n" + body

    out: list[str] = []
    paragraph: list[str] = []
    in_math = False
    math_lines: list[str] = []
    in_enum = False
    figure_index = 0

    def flush_paragraph() -> None:
        if paragraph:
            out.append(escape_text(" ".join(x.strip() for x in paragraph)))
            out.append("")
            paragraph.clear()

    def close_enum() -> None:
        nonlocal in_enum
        if in_enum:
            out.append(r"\end{enumerate}")
            out.append("")
            in_enum = False

    for raw in body.splitlines():
        line = raw.rstrip()
        if line.strip() == "$$":
            flush_paragraph(); close_enum()
            if in_math:
                out.extend([r"\begin{equation}", "\n".join(math_lines), r"\end{equation}", ""])
                math_lines.clear()
            in_math = not in_math
            continue
        if in_math:
            math_lines.append(line)
            continue
        heading = re.match(r"^(##|###)\s+(.+)$", line)
        if heading:
            flush_paragraph(); close_enum()
            title = heading.group(2).strip()
            command = "section" if heading.group(1) == "##" else "subsection"
            out.append(rf"\{command}{{{escape_text(title)}}}")
            out.append("")
            if title in FIGURES:
                figure_index += 1
                out.append(figure_tex(*FIGURES[title], figure_index))
            continue
        item = re.match(r"^\d+\.\s+(.+)$", line)
        if item:
            flush_paragraph()
            if not in_enum:
                out.append(r"\begin{enumerate}")
                in_enum = True
            out.append(r"\item " + escape_text(item.group(1)))
            continue
        if not line.strip():
            flush_paragraph(); close_enum()
        else:
            paragraph.append(line)
    flush_paragraph(); close_enum()
    return abstract, keywords, "\n".join(out)


def main() -> None:
    abstract, keywords, body = convert(SOURCE.read_text(encoding="utf-8"))
    document = PREAMBLE + rf"""
\begin{{abstract}}
{abstract}
\end{{abstract}}
\begin{{IEEEkeywords}}
{escape_text(keywords)}
\end{{IEEEkeywords}}
{body}
\end{{document}}
"""
    TEX.write_text(document, encoding="utf-8")
    BUILD.mkdir(exist_ok=True)
    subprocess.run(
        [str(TECTONIC), "-X", "compile", str(TEX), "--outdir", str(BUILD)],
        cwd=PROJECT,
        env={"PATH": f"{TECTONIC.parent}:{os.environ.get('PATH', '')}",
             "XDG_CACHE_HOME": "/tmp/scvicar-tectonic-cache"},
        check=True,
    )
    built = BUILD / TEX.with_suffix(".pdf").name
    VECTOR_PDF.write_bytes(built.read_bytes())
    source_pdf = fitz.open(built)
    compatible = fitz.open()
    scale = 160 / 72
    for page in source_pdf:
        pixmap = page.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False)
        image = pixmap.tobytes("jpeg", jpg_quality=90)
        target = compatible.new_page(width=page.rect.width, height=page.rect.height)
        target.insert_image(target.rect, stream=image)
    compatible.save(PDF, deflate=True, garbage=4)
    ASCII_PDF.write_bytes(PDF.read_bytes())
    print(PDF)


if __name__ == "__main__":
    main()
