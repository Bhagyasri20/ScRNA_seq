# Stage-wise metabolic modelling of gastric lesions

Build one genome-scale metabolic model for each stage of gastric carcinogenesis (NAG → CAG → IM → EGC) from single-cell RNA-seq, with nutrient availability set by Human Plasma-Like Medium (HPLM), and compare the models to find reactions and pathways that change between stages.

> **Status:** the pipeline is implemented. It has not yet been run end to end on the full data, so no results are included. Update this section when you add results.

## Overview

```
GSE134520 scRNA-seq (13 biopsies)
   │  Seurat: QC, doublet removal, Harmony, annotation, epithelial selection
   ▼
Epithelial pseudobulk per sample → CPM → averaged per stage
   │
   ▼
Recon3D+ (human metabolic model) + HPLM uptake bounds + gene-protein-reaction rules
   │  core reactions and weights from expression
   ▼
SPECTRA (spectraCCME, minNetLP) → one model per stage
   │
   ▼
FVA → flux span ratio (FSR) → subsystem enrichment
```

| Stage | Samples |
|---|---|
| NAG | NAG1–3 |
| CAG | CAG1–3 |
| IM | IMW1–2 (mild), IMS1–4 (severe), pooled |
| EGC | EGC (one sample) |

## Repository contents

| File | Language | Purpose |
|---|---|---|
| `Seurat_final_modified.R` | R | Single-cell processing from raw counts to annotated epithelial cells, and export of the stage-level inputs |
| `make_model_inputs.R` | R | Re-creates the MATLAB input files from the saved object `seu_annotated_v3.rds` (use if you only need the export step) |
| `stage_models.m` | MATLAB | HPLM bounds, SPECTRA extraction of the stage models, FVA, FSR and pathway enrichment |
| `make_hplm_templates.m` | MATLAB | Writes blank HPLM mapping/concentration templates and checks exchange IDs. **It overwrites `hplm_*.csv`: do not rerun it once your files are filled** |
| `build_mapping.py` | Python | Builds `hplm_to_exchange.csv` and `hplm_mapping_notes.csv` from the concentration file |
| `hplm_concentrations.csv` | data | HPLM concentrations (mmol/L) |
| `hplm_to_exchange.csv` | data | HPLM component → Recon3D+ exchange reaction |
| `hplm_mapping_notes.csv` | data | How each HPLM component was handled (mapped, salt form, ions, unmapped) |

The scripts read their inputs from the working folder, and the paths are set at the top of each file. The simplest layout is one flat folder containing all of the files above plus the downloaded data.

## Requirements

- **R** with Seurat (v5), harmony, scDblFinder, SingleCellExperiment, data.table, Matrix, dplyr, ggplot2, AnnotationDbi and org.Hs.eg.db
- **MATLAB** with the [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/) and a supported LP/MILP solver (the SPECTRA README lists CPLEX or Gurobi)
- **SPECTRA** code: <https://github.com/bisect-group/SPECTRA>
- **Python 3** (only for `build_mapping.py`)

## Data (not included)

| Item | Source |
|---|---|
| scRNA-seq counts | GEO accession **GSE134520** (Zhang et al., 2019); place `GSE134520_RAW.tar` in the working folder |
| Recon3D+ model | The SPECTRA repository (`Recon3D+/UpdatedRecon3D.mat`) or the SPECTRA paper's supplementary data |
| Recon3D+ gene table | `recon3dplus_gene_symbols_correct.csv` with columns `recon_id`, `symbol`, `entrez_id` |
| HPLM composition | Cantor et al. 2017, *Cell* 169:258 (supplementary composition table) |

## How to run

1. **Single-cell processing (R).** Set `base_dir` at the top of `Seurat_final_modified.R` and run it from top to bottom. It writes `model_input/N_per_stage.csv`, `stage_cpm_for_matlab.csv`, `sample_cpm_for_matlab.csv` and `sample_stage_map.csv`. If `seu_annotated_v3.rds` already exists, run `make_model_inputs.R` alone.
2. **HPLM files.** Use the provided `hplm_concentrations.csv` and `hplm_to_exchange.csv`. Check that the exchange IDs exist in your model (`stage_models.m` matches IDs regardless of punctuation style, for example `EX_glc_D(e)` and `EX_glc__D_e`). Rows with no concentration or no exchange ID are skipped.
3. **MATLAB setup.** In `stage_models.m`, set `cfg.spectraDir` to the SPECTRA folder, check that `which spectraCCME` finds the function, and set a solver with `changeCobraSolver`.
4. **Run** `stage_models` from the working folder. Output goes to `stage_results/<run>/`.

## Method summary

- **Expression profiles.** Counts of all epithelial cells are summed per sample, converted to CPM and averaged within each stage, so each biopsy has equal weight.
- **Medium.** All uptake is closed, then HPLM components are opened with
  `max uptake (mmol/gDW/h) = [C (mmol/L) × V (mL) / 1000] / [t (h) × DW (gDW/cell) × N (cells)]`,
  applied as `lb = −min(cap, value)`. Inorganic ions and gases are kept open.
- **Core reactions.** Called per sample, then kept for a stage if core in at least half of that stage's samples. The default rule is a fallback (top quartile of expressed reactions); the SPECTRA authors' LocalGini method is in the SPECTRA repository and can replace it.
- **Extraction.** `spectraCCME` with `minNetLP` and stoichiometric consistency; weights are 1/(1+normalised expression) for non-core reactions.
- **Comparison.** FVA (no optimality requirement), FSR with the thresholds of the SPECTRA paper (up: higher maximum flux and FSR > 2; down: lower maximum flux and FSR < 0.5), and hypergeometric subsystem enrichment with Benjamini–Hochberg correction.

## Parameters and assumptions

The biopsy dataset has no culture parameters, so these are **modelling assumptions** shared by all stages:

| Parameter | Value | Reason |
|---|---|---|
| V | 10 mL | common dish/flask volume |
| t | 48 h | usual medium-change interval |
| DW | 3×10⁻¹⁰ gDW/cell | mid-range mammalian cell dry weight (verify against literature) |
| N | 10⁶ cells | culture-scale cell number; gives glucose ≈ 3.5 mmol/gDW/h, a plausible upper bound |
| cap | 1000 | prevents unrealistic bounds for abundant ions |

Using the number of sequenced cells as N would loosen the bounds roughly 25–100 times, so it is run only as a sensitivity analysis. Four runs are performed: the main run, per-stage N, and the medium scaled by 0.5 and by 2.

## Outputs

`stage_results/<run>/` contains `model_<stage>.mat`, `fva_<stage>.csv`, `blocked_core_<stage>.csv`, `up_*.csv`, `down_*.csv` and `enrichment_*.csv` for every comparison, and `stage_results/summary_all_runs.csv` with counts across runs.

## Limitations

- HPLM models plasma, whereas gastric mucosa is supplied by tissue fluid.
- V, t, DW and N are assumptions; uptake bounds are upper limits, not measured fluxes.
- EGC is a single sample, so it has no replicate support; mild IM has only two samples.
- IM pools mild and severe samples and is therefore weighted towards severe IM; H. pylori status is partly confounded with severity (IMW1, IMS1 and IMS2 are positive).
- Harmony correction by sample can remove real biology for a stage that is a single sample.
- The default core-reaction rule is a fallback, and no experimental validation or comparison with another method is included.

## References

1. Kumar S P, et al. *Generalist method to reconstruct metabolic networks from multi-omics data at large-scale* (SPECTRA). bioRxiv 2026. doi:10.64898/2026.04.02.716249
2. Zhang P, et al. Cell Reports 2019;27:1934 (GSE134520).
3. Cantor JR, et al. Cell 2017;169:258–272 (HPLM).
4. Brunk E, et al. Nature Biotechnology 2018;36:272–281 (Recon3D).
5. Heirendt L, et al. Nature Protocols 2019;14:639–702 (COBRA Toolbox).
6. Kumar S P, Bhatt NP. npj Systems Biology and Applications 2025;11:142 (LocalGini).
7. Korsunsky I, et al. Nature Methods 2019;16:1289 (Harmony).

## Contact

Contact: bt25m016@smail.iitm.ac.in
