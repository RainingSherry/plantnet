from __future__ import annotations

import hashlib
import html
import re
from pathlib import Path

import markdown
from matplotlib import mathtext
from weasyprint import HTML


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "scVICAR_中文阅读稿.md"
BUILD = ROOT / "build"
MATH_DIR = BUILD / "math"
HTML_OUT = BUILD / "scVICAR_中文阅读稿.html"
PDF_OUT = ROOT / "scVICAR_中文阅读稿.pdf"
FONT_PATH = Path.home() / ".fonts" / "NotoSansSC.ttf"


CSS = r"""
@font-face {
  font-family: "Noto Sans SC";
  src: url("__FONT_URI__");
  font-weight: 100 900;
}
@page {
  size: A4;
  margin: 16mm 14mm 17mm 14mm;
  @top-left {
    content: "scVICAR 中文阅读稿";
    font: 7.5pt "Noto Sans SC";
    color: #555;
  }
  @top-right {
    content: "TCBB 风格双栏版";
    font: 7.5pt "Noto Sans SC";
    color: #555;
  }
  @bottom-center {
    content: counter(page);
    font: 8pt "Noto Sans SC";
  }
}
* { box-sizing: border-box; }
body {
  margin: 0;
  color: #111;
  font-family: "Noto Sans SC", sans-serif;
  font-size: 8.65pt;
  font-weight: 400;
  line-height: 1.48;
  text-align: justify;
}
.title-block {
  text-align: center;
  margin: 3mm 8mm 5mm;
}
h1 {
  margin: 0 0 3mm;
  font-size: 18pt;
  line-height: 1.25;
  font-weight: 700;
  letter-spacing: .1pt;
}
.authors { font-size: 9pt; margin-bottom: 4mm; }
.abstract {
  border-top: .55pt solid #222;
  border-bottom: .55pt solid #222;
  padding: 2.5mm 4mm 2.7mm;
  margin: 0 5mm 5mm;
  font-size: 8.4pt;
  line-height: 1.46;
  text-align: justify;
}
.abstract h2 {
  display: inline;
  margin: 0;
  font-size: 8.6pt;
  font-weight: 700;
}
.abstract h2::after { content: "—"; }
.abstract p { display: inline; margin: 0; }
.abstract p + p { display: block; margin-top: 1.3mm; }
.paper-body {
  column-count: 2;
  column-gap: 6.5mm;
  column-rule: .25pt solid #ddd;
}
h2, h3, h4 { break-after: avoid; page-break-after: avoid; }
h2 {
  margin: 3.1mm 0 1.5mm;
  font-size: 10.2pt;
  line-height: 1.25;
  font-weight: 700;
}
h3 {
  margin: 2.5mm 0 1.1mm;
  font-size: 9.2pt;
  line-height: 1.25;
  font-weight: 700;
}
h4 {
  margin: 2mm 0 1mm;
  font-size: 8.8pt;
  font-weight: 700;
}
p { margin: 0 0 1.55mm; orphans: 3; widows: 3; }
ol, ul { margin: 1mm 0 1.8mm 4.5mm; padding-left: 3.4mm; }
li { margin: 0 0 .8mm; }
strong { font-weight: 700; }
code {
  font-family: "DejaVu Sans Mono", monospace;
  font-size: .88em;
  background: #f3f3f3;
  padding: 0 .5mm;
}
.equation {
  text-align: center;
  margin: 1.8mm 0 2mm;
  break-inside: avoid;
}
.equation img {
  display: inline-block;
  max-width: 94%;
  max-height: 21mm;
}
.inline-math {
  display: inline-block;
  height: 1.15em;
  width: auto;
  vertical-align: -0.22em;
  margin: 0 .12em;
}
table {
  border-collapse: collapse;
  width: 100%;
  margin: 2mm 0;
  font-size: 7.8pt;
  break-inside: avoid;
}
th, td { border: .4pt solid #777; padding: .8mm 1mm; }
th { font-weight: 700; background: #f2f2f2; }
.footer-note {
  margin-top: 4mm;
  padding-top: 2mm;
  border-top: .4pt solid #777;
  font-size: 7.5pt;
  color: #555;
}
"""


def render_math(expression: str, display: bool) -> Path:
    normalized = expression.strip().replace("\n", " ")
    normalized = re.sub(r"\\mathcal\s+([A-Za-z])", r"\\mathcal{\1}", normalized)
    normalized = re.sub(r"\\mathbb\s+([A-Za-z])", r"\\mathbb{\1}", normalized)
    normalized = re.sub(r"\\mathbf\s+([A-Za-z0-9])", r"\\mathbf{\1}", normalized)
    normalized = re.sub(r"\\(?:big|Big|bigg|Bigg)[lr]?", "", normalized)
    normalized = re.sub(r"\\le(?![A-Za-z])", r"\\leq", normalized)
    normalized = re.sub(r"\\ge(?![A-Za-z])", r"\\geq", normalized)
    digest = hashlib.sha256((("D" if display else "I") + normalized).encode()).hexdigest()[:16]
    output = MATH_DIR / f"{digest}.svg"
    if output.exists():
        return output
    wrapped = f"${normalized}$"
    try:
        mathtext.math_to_image(wrapped, output, format="svg", dpi=220)
    except Exception:
        fallback = re.sub(r"\\(?:left|right)", "", normalized)
        fallback = fallback.replace(r"\operatorname", r"\mathrm")
        mathtext.math_to_image(f"${fallback}$", output, format="svg", dpi=220)
    return output


def protect_math(text: str) -> tuple[str, dict[str, str]]:
    replacements: dict[str, str] = {}

    def display(match: re.Match[str]) -> str:
        token = f"MATHDISPLAY{len(replacements):04d}TOKEN"
        path = render_math(match.group(1), True)
        replacements[token] = (
            f'<div class="equation"><img src="{path.as_uri()}" alt="'
            f'{html.escape(match.group(1).strip())}"></div>'
        )
        return f"\n\n{token}\n\n"

    text = re.sub(r"\$\$(.*?)\$\$", display, text, flags=re.S)

    def inline(match: re.Match[str]) -> str:
        token = f"MATHINLINE{len(replacements):04d}TOKEN"
        path = render_math(match.group(1), False)
        replacements[token] = (
            f'<img class="inline-math" src="{path.as_uri()}" alt="'
            f'{html.escape(match.group(1).strip())}">'
        )
        return token

    text = re.sub(r"(?<!\$)\$([^\n$]+?)\$(?!\$)", inline, text)
    return text, replacements


def markdown_html(text: str) -> str:
    protected, replacements = protect_math(text)
    rendered = markdown.markdown(
        protected,
        extensions=["extra", "sane_lists"],
        output_format="html5",
    )
    for token, replacement in replacements.items():
        rendered = rendered.replace(f"<p>{token}</p>", replacement)
        rendered = rendered.replace(token, replacement)
    return rendered


def main() -> None:
    BUILD.mkdir(parents=True, exist_ok=True)
    MATH_DIR.mkdir(parents=True, exist_ok=True)
    source = SOURCE.read_text(encoding="utf-8")
    source = re.sub(r"^# .*?\n", "", source, count=1)
    source = re.sub(r"^> .*?\n", "", source, count=1, flags=re.M)
    source = source.split("\n## 作者修改提示", 1)[0].strip()
    abstract_marker = "## 摘要"
    intro_marker = "## 1. 引言"
    abstract_start = source.index(abstract_marker)
    intro_start = source.index(intro_marker)
    abstract_md = source[abstract_start:intro_start].strip()
    body_md = source[intro_start:].strip()

    css = CSS.replace("__FONT_URI__", FONT_PATH.as_uri())
    document = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<title>scVICAR 中文阅读稿</title><style>{css}</style></head><body>
<header class="title-block">
  <h1>scVICAR：用于掩码单细胞表征学习的有界图邻域锚点恢复</h1>
  <div class="authors">匿名作者 · 双盲评审中文阅读版</div>
</header>
<section class="abstract">{markdown_html(abstract_md)}</section>
<main class="paper-body">{markdown_html(body_md)}
<p class="footer-note">本中文稿用于作者阅读与论证核对；投稿与引用以英文正式稿为准。</p>
</main></body></html>"""
    HTML_OUT.write_text(document, encoding="utf-8")
    HTML(filename=str(HTML_OUT), base_url=str(ROOT)).write_pdf(PDF_OUT)
    print(PDF_OUT)


if __name__ == "__main__":
    main()
