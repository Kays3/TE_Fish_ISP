# Manifests

| File | What it lists | How it was made |
|---|---|---|
| `inputs.sha256.tsv` | The large inputs that stay outside git: 162 TEcount tables of the *A. polyacanthus* cohort, the annotation, and the zebrafish counts, DESeq2 tables and alignment QC. Columns: sha256, bytes, path below `DATA_ROOT`. | Rows copied from the verified transfer manifest of the archive that now holds these folders (sha256 recorded at copy time, 2026-10-05; manifest sha256 `0d1022628e2796182482e8fc50b3d31b6c5eb97e63150fe0fb80e82cf34a398d`). No file was read or copied through a laptop. |
| `repo_files.sha256` | Every file of this repository that other files depend on: `data/small_inputs/`, `config/`, `results/`. | `scripts/write_checksums.sh` |
| `figures.sha256` | The committed figures (PNG). | `scripts/write_checksums.sh` |

`DATA_ROOT` is the directory that contains `te-oa-full/`. The raw reads (834.6 GB; accessions in
`config/zebrafish_perturbation_studies.tsv`, and every FASTQ with its ENA MD5 in
`pipelines/zebrafish_validation/zf_download_list.tsv`) are not part of the manifest; they are re-downloadable and
exist in one archived copy.

Check with `scripts/verify_inputs.sh` (repository files always; `DATA_ROOT` files when the variable is set).
