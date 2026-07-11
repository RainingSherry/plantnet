# scVICAR 中文阅读稿

本目录用于作者内部阅读、核对和修改，不作为投稿文件。

- `scVICAR_中文阅读稿.md`：与英文主文结构对应的完整中文稿。
- `scVICAR_中文阅读稿.tex`：由 Markdown 生成的 IEEEtran 中文 LaTeX 稿。
- `scVICAR_中文阅读稿.pdf`：无字体依赖的兼容交付版，复用英文主文六幅图，可在不同 PDF 查看器中稳定显示中文。
- `scVICAR_zh_reading.pdf`：与兼容交付版内容相同的纯 ASCII 文件名副本，用于不支持中文路径的预览器。
- `scVICAR_中文阅读稿_vector.pdf`：Tectonic 直接生成的可搜索矢量版，供构建归档使用；部分查看器可能不兼容其 Fandol CID 字体映射。
- `build_chinese_latex.py`：生成 LaTeX 并编译 PDF 的可重复构建脚本。
- 英文稿仍是唯一正式版本；中文修改确认后，应同步回英文源文件。
- 公式、实验数值、统计限定和结论边界均按英文稿保留。

## 构建

在项目根目录运行：

```bash
/data/luolie/conda/base/bin/python papers/scVICAR/manuscript_zh/build_chinese_latex.py
```

脚本固定调用仓库内 `.codex_tex/bin/tectonic`，并将 Tectonic 缓存写入
`/tmp/scvicar-tectonic-cache`。`build_latex/` 保存编译产物，目录根部的
`scVICAR_中文阅读稿.pdf` 和 `scVICAR_zh_reading.pdf` 是跨查看器交付版本。脚本同时保留
`scVICAR_中文阅读稿_vector.pdf`，然后以 160 dpi 将页面固化为不依赖字体映射的兼容 PDF。

## 术语约定

| 英文术语 | 中文阅读稿用法 |
|---|---|
| graph-vicinal anchor recovery | 图邻域锚点恢复 |
| topology-informed affinity | 拓扑信息亲和度 |
| cell-specific mixing coefficient | 逐细胞混合系数 |
| matched-backbone | 骨干匹配 / matched-backbone |
| frozen protocol | 冻结协议 |
| NoMix | NoMix，不混合对照 |
| marker coherence | 标记基因一致性 |
| low-label transductive annotation | 低标签量传导式注释 |
