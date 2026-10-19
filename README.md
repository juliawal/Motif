# Motif

**Motif** (**M**ethylation fact**o**r iden**tif**ication) identifies candidate methylation factors from paired TCGA gene expression and DNA methylation data (UCSC Xena, GDC Pan-Cancer). The pipeline combines GRNBoost2-based network inference, filtering, stability selection, and evaluation of the selected candidates.

## Installation

The pipeline uses Conda and Slurm on Linux. Run all commands from the project
directory.

```bash
conda env create -f environment.yml
conda activate motif
```

Set the cohorts, Slurm partition and account in [config.yaml](config.yaml).
The paths below match this configuration.

## Input data

Download these four files from the UCSC Xena GDC Pan-Cancer dataset to `data/raw/`:

| File | Contents |
| --- | --- |
| [GDC-PANCAN.htseq_fpkm-uq.tsv](https://gdc.xenahubs.net/download/GDC-PANCAN.htseq_fpkm-uq.tsv) | Gene expression; log2(FPKM-UQ + 1) |
| [GDC-PANCAN.htseq_counts.tsv](https://gdc.xenahubs.net/download/GDC-PANCAN.htseq_counts.tsv) | Counts stored as log2(count + 1) |
| [GDC-PANCAN.methylation450.tsv](https://gdc.xenahubs.net/download/GDC-PANCAN.methylation450.tsv) | 450K methylation beta values |
| [GDC-PANCAN.basic_phenotype.tsv](https://gdc.xenahubs.net/download/GDC-PANCAN.basic_phenotype.tsv) | Sample IDs, project, sample type, sex and age |

Download all four files with:

```bash
mkdir -p data/raw resources
for name in htseq_fpkm-uq htseq_counts methylation450 basic_phenotype; do
  curl -fL "https://gdc.xenahubs.net/download/GDC-PANCAN.${name}.tsv" \
    -o "data/raw/GDC-PANCAN.${name}.tsv" || break
done
```

Run the preparation steps in this order:

```bash
# Keep samples shared by all three matrices and usable metadata; split by TCGA project.
python -m pipeline.split_data_to_cohorts \
  --metadata data/raw/GDC-PANCAN.basic_phenotype.tsv \
  --expression data/raw/GDC-PANCAN.htseq_fpkm-uq.tsv \
  --methylation data/raw/GDC-PANCAN.methylation450.tsv \
  --raw-counts data/raw/GDC-PANCAN.htseq_counts.tsv \
  --output-dir splitted_by_project

# Move files into data/<cohort>/ for the cohorts in config.yaml.
python -m scripts.fill_cohort_dirs --config config.yaml \
  --split-dir splitted_by_project --data-dir data

# Convert log2(count + 1) back to rounded integer counts; overwrite these files.
python -m scripts.repair_raw_counts data/*/raw_counts.tsv

# Save sex metadata for C1 and age/sex metadata for F4.
python -m scripts.build_metadata \
  --metadata data/raw/GDC-PANCAN.basic_phenotype.tsv \
  --config config.yaml --data-directory data
```

Each `data/<cohort>/` now contains `expression.tsv`, `raw_counts.tsv`,
`methylation.tsv`, `shared_samples.tsv`, `metadata.tsv` and `age_and_sex.tsv`.
SKCM files are copied into `SKCM_PRIMARY` and `SKCM_METASTATIC`; stage A selects
the sample types and removes patients present in both groups.

## Reference files

Download the linked files and save them in `resources/` under the names in the
first column. Keep the `.gz` files compressed except for the genome, which is
unpacked below. Data and references are not included in this repository.

| Save as in `resources/` | What it is and where to download it |
| --- | --- |
| `gencode.v22.annotation.gtf.gz` | Gene and transcript positions: [GENCODE v22 GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz) |
| `hgnc_complete_set.txt` | Gene symbols, aliases and Entrez IDs: [HGNC complete set](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt) |
| `HM450.hg38.manifest.gencode.v36.probeMap` | 450K probe coordinates: [UCSC Xena probeMap](https://gdc.xenahubs.net/download/HM450.hg38.manifest.gencode.v36.probeMap) |
| `HM450.hg38.manifest.zhou.tsv.gz` | Probe coordinates and `MASK_general` for removing masked probes: [Zhou HM450 manifest](https://zhouserver.research.chop.edu/InfiniumAnnotation/current/HM450/HM450.hg38.manifest.tsv.gz) |
| `HM450.hg38.manifest.gencode.v36.zhou.tsv.gz` | CpG island context for the F1 annotation builder: [Zhou GENCODE v36 manifest](https://zhouserver.research.chop.edu/InfiniumAnnotation/current/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz) |
| `hg38.fa.gz` | Human genome sequence: [UCSC hg38 FASTA](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) |
| `JASPAR2024_CORE_vertebrates.meme` | DNA-binding motifs for FIMO: [JASPAR2024 CORE vertebrates, non-redundant MEME](https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt) |
| `JASPAR2024.sqlite` | Gene-to-motif lookup: [JASPAR2024 database from pyjaspar 4.0.0](https://raw.githubusercontent.com/asntech/pyjaspar/v4.0.0/pyjaspar/data/JASPAR2024.sqlite) |
| `hippie.txt` | Human protein interactions: [HIPPIE v2.3 TAB](https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_v2_3.txt) |

Build the reference tables and download the gene-set libraries:

```bash
# Build gene information, transcript TSS positions and the protein-coding gene list.
python -m scripts.build_gene_annotation \
  --gencode-file resources/gencode.v22.annotation.gtf.gz \
  --hgnc-file resources/hgnc_complete_set.txt \
  --gene-annotation-output resources/gene_annotation.v22.tsv \
  --tss resources/gene_annot.tsv \
  --coding-genes resources/gencode.v22.protein_coding_gene_ids.txt

# Build CpG coordinates for C2 and F4.
python -m scripts.build_cpg_coordinates \
  --cpg-map resources/HM450.hg38.manifest.gencode.v36.probeMap \
  --output-file resources/cpg_annot.tsv

# Keep CpG island context; calculate gene overlaps and TSS distances with GENCODE v22.
python -m scripts.build_cpg_annotations_for_f1 \
  --cpg-annotations-file resources/HM450.hg38.manifest.gencode.v36.zhou.tsv.gz \
  --gencode-file resources/gencode.v22.annotation.gtf.gz \
  --output resources/cpg_gene_annotations.hg38.gencode.v22.tsv.gz

# Unpack the genome and prepare it for reading DNA sequences around CpGs.
gunzip resources/hg38.fa.gz
python -m scripts.build_fasta_index resources/hg38.fa

# Download the libraries listed in config.yaml as GMT files.
python -m scripts.download_gene_set_enrichment_libraries \
  --config config.yaml --output-dir resources/gene_sets
```

The libraries come from [Enrichr](https://maayanlab.cloud/Enrichr/):
`MSigDB_Hallmark_2020`, `Reactome_2022`, `GO_Biological_Process_2023` and
`GO_Molecular_Function_2023`. ChIP-seq peaks are downloaded after F0, when the
selected MFs are known (see below).

## Run

Run the complete pipeline on Slurm:

```bash
motif all
```

This runs A to F0, downloads [ReMap 2022](https://remap.univ-amu.fr/download_page)
ChIP-seq peaks for selected MFs without a JASPAR motif, then runs F1 to G.
If a run stops, call `motif all` again with the same inputs and configuration.

To run one stage, use its name, for example `motif A`, `motif D` or `motif F2`.
These commands run only that stage; its input files from earlier stages must
already exist. To run everything one stage at a time:

```bash
motif A
motif B
motif C1
motif C2
motif D
motif E
motif F0
python -m scripts.download_chipseq_peaks --from-pipeline-results results --config config.yaml
motif F1
motif F2
motif F3
motif F4
motif G
```

For execution on the current machine instead of Slurm, add `--executor local`
and a core limit, for example `motif all --executor local --cores 32` or
`motif D --executor local --cores 4`.

## Stages

| Stage | What it does |
| --- | --- |
| A | Select tumour samples, filter missing data and masked CpGs, select highly variable genes and build the feature matrix. |
| B | Run GRNBoost2 with several seeds to score gene–CpG edges. |
| C1 | Remove edges where both gene expression and CpG methylation depend on sex. |
| C2 | Remove negatively correlated edges where the CpG is near the MF gene's TSS (reverse-causality filter). |
| D | Compare MF lists from split halves of the GRN runs using Kuncheva stability; rank settings across cohorts. |
| E | Apply each selected setting to all runs, aggregate edges and select the top MFs per cohort. |
| F0 | Add gene symbols and prepare the selected MF–CpG pairs for evaluation. |
| F1 | Compare motif hits or ChIP-seq overlaps at target CpGs with matched background CpGs. |
| F2 | Measure proximity between selected MFs in the HIPPIE protein-protein-interaction network. |
| F3 | Test gene-set enrichment of the selected MFs. |
| F4 | Test MF expression to CpG methylation to nearby gene expression using mediation models with age and sex metadata. |
| G | Collect evaluation results across all cohorts for each selected setting. |

## Results

**A to D: per cohort**, under `results/<cohort>/`:

| Stage | Files |
| --- | --- |
| A | `prepared/feature_matrix.tsv`, `regulator_names.tsv`, `target_names.tsv`, `samples.tsv`, `counts.tsv` and `cpg_region_map.tsv` (all in `prepared/`) |
| B | `grn/seed_<seed>.tsv`: edges and weights for each run |
| C1 | `filtered/after_sex.tsv`, `gene_sex_test.tsv` and `cpg_sex_test.tsv` (all in `filtered/`) |
| C2 | `filtered/retained_edges.tsv`: pairs retained after both filters |
| D | `stability/grid_stability.tsv`: stability of every tested setting |

**D: comparison across cohorts**, under `results/stability_comparison/`:

| File | Contents |
| --- | --- |
| `ranking.tsv` | Settings ranked by median Kuncheva stability across cohorts |
| `selection.tsv` | The selected best settings and their ranks |
| `kuncheva_by_cohort.tsv` | Stability matrix: cohorts x settings |
| `selected_stability.tsv` | Stability per cohort for the selected settings |
| `stability_results_with_enough_mfs.tsv` | Settings that produced enough MFs in every cohort and split half |

A setting specifies the edge cut, run aggregation and MF scoring method.
`rank_1` is the setting with the highest median stability across cohorts.
The current configuration keeps three settings (`n_selected_stable_settings`).
Each top setting is used for every cohort.

**E to F4: per setting rank**, under `results/selected_settings/rank_<rank>/`:

| Stage | Path within the rank folder |
| --- | --- |
| E | `selection_applied.tsv`: the setting used for this rank |
| E | `<cohort>/aggregated/edges.tsv`: aggregated edges |
| E | `<cohort>/selected_mfs/mfs.tsv`: selected MFs and scores |
| F0 | `<cohort>/evaluation_input/mf_cpg_pairs.tsv` and `mf_mapping.tsv` |
| F1 | `<cohort>/evaluation/binding_sites.tsv`, `binding_sites_per_mf.tsv` and `binding_evidence/` with FIMO output and CpG assignments |
| F2 | `<cohort>/evaluation/ppi_proximity.tsv` |
| F3 | `<cohort>/evaluation/gene_set_enrichment.tsv` |
| F4 | `<cohort>/evaluation/mediation.tsv` and `mediation_triplets.tsv` |

**G: per setting rank, across cohorts**, under
`results/selected_settings/rank_<rank>/evaluation_summary/`:

- `evaluation_by_cohort.tsv`: combined evaluation results for each cohort.
- `evaluation_summary.tsv`: counts of testable cohorts and full MF lists for each evaluation stage.
- `enrichment_top_terms.tsv`: the top enrichment terms per cohort, selected across libraries.

Use the saved results to make plots and manuscript tables separately by hand.

The GRNBoost2 code in [adapted_grnboost/](adapted_grnboost/) is adapted from Arboreto.
The package lives directly in the project directory and is imported as
`adapted_grnboost`. Its license is kept in
[LICENSE.upstream.txt](adapted_grnboost/LICENSE.upstream.txt).
