# Data provenance and license audit

Updated: 2026-07-10. Status: **closed** for the private analytical reuse and
publication described here. This does not authorize public redistribution of
the source H5AD files.

The local H5AD files do not contain license or citation fields in `.uns`.
Therefore public availability is not treated as permission, and no license is
inferred from dataset names alone.

| scVICAR dataset | Local provenance identifier | Authoritative license evidence | Status |
|---|---|---|---|
| Blood_BoneMarrow | CELLxGENE `dataset_id=d3566d6a-a455-4a15-980f-45eb29114cab`; Triana et al., DOI `10.1038/s41590-021-01059-0` | Census 2024-07-01 citation resolves the collection; Europe PMC full text states CC BY 4.0 | Resolved |
| Human_Pancreas_1 | CELLxGENE `dataset_id=66d15835-5dc8-4e96-b0eb-f48971cb65e8`; Enge et al., DOI `10.1016/j.cell.2017.09.004`; GEO `GSE81547` | Census citation and primary data-availability statement verified; NCBI places no restrictions on molecular-data use/distribution but cannot transfer any submitter IP rights | Resolved for analysis; no public redistribution claim |
| Human_Pancreas_3 | Baron et al., DOI `10.1016/j.cels.2016.08.011`; GEO `GSE84133` | Primary data-availability statement and NCBI molecular-data policy verified | Resolved for analysis; no public redistribution claim |
| Mouse_Pancreas_1 | Baron et al., DOI `10.1016/j.cels.2016.08.011`; GEO `GSE84133` | Primary data-availability statement and NCBI molecular-data policy verified | Resolved for analysis; no public redistribution claim |
| PRJNA895163 | NCBI BioProject `PRJNA895163`; Song et al., *Frontiers in Plant Science* 2022, DOI `10.3389/fpls.2022.1053669` | Europe PMC full text states Creative Commons Attribution (CC BY); attribution and article citation required | Resolved |
| TabulaSapiens_Pancreas | CELLxGENE `dataset_id=ff45e623-7f5f-46e3-b47d-56be0341f66b`; Tabula Sapiens count/metadata release DOI `10.6084/m9.figshare.14267219.v5` | Official Figshare API reports CC BY 4.0; article distinguishes public count/metadata from raw reads requiring a data-transfer agreement | Resolved for the count/metadata derivative used here |

## Checks performed

- Inspected `.uns` and relevant `.obs` provenance fields in all six source
  H5ADs.
- Recorded the three available CELLxGENE dataset UUIDs and the BioProject
  accession.
- Attempted a CELLxGENE Census metadata query against release `2025-11-08`;
  the dataset table exposes citation/collection DOI fields but not a license
  column, and the remote S3 query failed checksum validation before records
  could be retrieved.
- Resolved `PRJNA895163` through the official NCBI BioProject summary and the
  originating open-access article. Europe PMC's archival full text contains an
  explicit CC BY license statement.
- Resolved the Tabula Sapiens count/metadata release through the primary
  article's data-availability statement and the official Figshare API. The
  analysis uses expression counts/metadata, not controlled raw sequence reads.
- Recovered the Blood/Bone-Marrow, Human Pancreas 1, and Tabula Sapiens source
  citations from the official CELLxGENE Census 2024-07-01 dataset table.
- Verified the Human Pancreas 1 and Baron human/mouse pancreas GEO accessions in
  their primary articles. NCBI's molecular-data policy states that NCBI places
  no restrictions on use or distribution, while warning that it cannot transfer
  possible submitter rights. Consequently scVICAR claims analytical reuse, not
  a right to republish the source H5ADs.

## Authoritative links resolved so far

- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA895163
- https://doi.org/10.3389/fpls.2022.1053669
- https://europepmc.org/articles/PMC9848496
- https://doi.org/10.1126/science.abl4896
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9812260/
- https://doi.org/10.6084/m9.figshare.14267219.v5
- https://api.figshare.com/v2/articles/14267219
- https://doi.org/10.1038/s41590-021-01059-0
- https://europepmc.org/articles/PMC8642243
- https://doi.org/10.1016/j.cell.2017.09.004
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81547
- https://doi.org/10.1016/j.cels.2016.08.011
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133
- https://www.ncbi.nlm.nih.gov/home/about/policies/#data

## Acceptance rule

Before submission, each row must link to an official collection, repository,
or primary-study page that states reuse terms. Any attribution, non-commercial,
controlled-access, or redistribution restriction must be reflected in the data
availability statement. The remote scVICAR archive remains private and must not
be advertised as a redistributable public dataset until this audit closes.
