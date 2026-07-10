# Name and citation audit

Checked: 2026-07-10 (Asia/Shanghai).

## scVICAR name screen

Exact-term searches returned no records for `scVICAR` in:

- Crossref Works (`query.title=scVICAR`, 0 results);
- Europe PMC (`query=scVICAR`, 0 results);
- arXiv (`all:scVICAR`, 0 results).

This is a publication-name screen, not a legal trademark clearance. The search
must be repeated immediately before submission because indexes can change.

## scCluBench citation

Crossref resolves the formal record as:

> Xu P, Wang Z, Wang Z, et al. scCluBench: Comprehensive Benchmarking of
> Clustering Algorithms for Single-Cell RNA Sequencing. Proceedings of the AAAI
> Conference on Artificial Intelligence. 2026;40(2):1364–1372.
> DOI: 10.1609/aaai.v40i2.37110.

The manuscript BibTeX entry matches the DOI, author list, year, volume, issue,
and pages. The paper abstract explicitly names marker-gene identification and
cell-type annotation as downstream biological tasks, supporting the limited
motivation cited in scVICAR. scVICAR does not claim to reproduce every
scCluBench task or its full benchmark ranking.

## Sources

- https://api.crossref.org/works?query.title=scVICAR&rows=5
- https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=scVICAR&format=json&pageSize=20
- https://export.arxiv.org/api/query?search_query=all%3AscVICAR&max_results=10
- https://api.crossref.org/works?query=scCluBench&rows=10
- https://doi.org/10.1609/aaai.v40i2.37110
