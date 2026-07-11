# scVICAR name-collision audit

Audit date: 2026-07-11 (Asia/Shanghai).

The exact token `scVICAR` and title/abstract variants were checked before the
paper snapshot was finalized.

| Source | Query | Result |
|---|---|---|
| Crossref REST API | title query `scVICAR` | 0 works |
| PubMed E-utilities | `scVICAR[Title/Abstract]` | 0 records |
| GitHub repository search API | repository query `scVICAR` | 0 repositories |
| PyPI JSON API | project `scvicar` | HTTP 404 (no project) |

OpenAlex fallback and the general web-search endpoint both failed at the
network layer and were not interpreted as negative results.  The successful
independent scholarly and software-index queries found no exact indexed
collision.  The name can therefore be retained for the draft, with a final
trademark and newly published-work check recommended immediately before public
release.  This audit does not claim legal clearance.
