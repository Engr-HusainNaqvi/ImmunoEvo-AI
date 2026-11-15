# ImmunoEvo-AI — Findings

Husain Naqvi · Project 14 · November 2026

Every number below comes from a CSV in `results/`, and the file is named next to the result.

## 1. Question and design

Immune-checkpoint blockade (ICB) response signatures are usually developed and tested inside one cohort. A clinically useful score has to work in a **new cohort, cancer type and drug**. I therefore evaluated everything with **leave-one-cohort-out (LOCO)** validation across 9 public ICB RNA-seq cohorts:
- melanoma ×4, renal ×2, urothelial, triple-negative breast, pancreatic;
- 820 pre-treatment patients, 256 responders.

**Expression.** Expression is z-scored within each cohort. No information from the held-out cohort, labels or otherwise, is used in training.

**Metric.**
- AUC is computed within each held-out cohort.
- The headline metric is the pooled AUC of within-cohort ranks. It weights patients rather than cohorts and avoids rewarding calibration differences between cohorts.
- 95% CIs come from a bootstrap stratified by cohort (resampling patients within each cohort).
- Paired comparisons use the same resamples for both models (B = 2,000; `loco_pairwise.csv`).

**Features.**
- **73 knowledge programmes:** 23 curated immune and TME sets (T-cell inflamed, CYT, IFN-γ, exhaustion, TLS/B cells, CAF/TGF-β, myeloid suppression, antigen presentation, etc.) plus the 50 MSigDB Hallmarks. Each is scored as the mean z-score of its genes.
- **Tumour evolution covariates** from the matched mutation calls:
  - TMB;
  - MATH (median-absolute-deviation VAF heterogeneity);
  - clonal fraction (cancer-cell fraction ≥ 0.75 after VAF-KDE purity correction);
  - clonal non-synonymous count;
  - each with a missingness flag.

**Models (all pre-specified before the LOCO results were seen):**
- **Fixed scores:** TIDE (the published binary responder call as distributed by iAtlas; the continuous score was not available, which handicaps its AUC), PD-L1/CD274, CD8A, CYT, GEP-18, IFNG-6, TMB.
- **Trained models:** logistic regression (L2, inner leave-one-cohort-out CV) on genes, on programmes, and on programmes + evolution covariates; random forest and gradient boosting on programmes.
- **ImmunoEvo-Net**, a knowledge-guided deep network:
  - a sparse gene→programme layer masked by programme membership;
  - a programme-level MLP;
  - evolution covariates concatenated;
  - a gradient-reversal cohort discriminator (λ = 0.3) that pushes the representation to be cohort-invariant;
  - a 10-seed ensemble.
- **Ablations:** no adversary, no evolution covariates, random programmes (same sizes, random genes).

## 2. Result 1: the simple knowledge model transfers best; the deep model does not transfer (`loco_summary.csv`)

| Model | Pooled AUC | 95% CI | Mean cohort AUC |
|---|---|---|---|
| LR programmes + evolution | **0.628** | 0.588–0.669 | 0.631 |
| TMB (n = 587, 8 cohorts) | 0.622 | 0.573–0.669 | 0.585 |
| IFNG-6 | 0.616 | 0.576–0.657 | 0.626 |
| CD8A | 0.615 | 0.573–0.654 | 0.641 |
| CYT | 0.610 | 0.567–0.651 | 0.612 |
| RF programmes | 0.610 | 0.569–0.650 | 0.599 |
| LR programmes | 0.609 | 0.569–0.652 | 0.610 |
| GEP-18 | 0.606 | 0.565–0.648 | 0.619 |
| GBM programmes | 0.590 | 0.548–0.631 | 0.565 |
| PD-L1 (CD274) | 0.589 | 0.549–0.630 | 0.591 |
| TIDE | 0.563 | 0.521–0.606 | 0.549 |
| LR genes (20,890) | 0.555 | 0.513–0.597 | 0.541 |
| ImmunoEvo-Net | 0.541 | 0.499–0.583 | 0.535 |
| Net w/o adversary | 0.540 | | |
| Net w/o evolution | 0.539 | | |
| Net, random programmes | 0.531 | | |

Paired differences for the programme + evolution model, with bootstrap 95% CI and two-sided p:

| Comparator | Δ AUC | 95% CI | p |
|---|---|---|---|
| ImmunoEvo-Net | +0.087 | 0.047, 0.125 | < 0.001 |
| LR genes | +0.073 | 0.041, 0.103 | < 0.001 |
| TIDE | +0.064 | 0.010, 0.122 | 0.023 |
| PD-L1 | +0.038 | 0.003, 0.073 | 0.029 |
| GBM programmes | +0.038 | 0.005, 0.070 | 0.027 |
| GEP-18 | +0.022 | 0.002, 0.042 | 0.027 |
| LR programmes (no evolution) | +0.019 | 0.007, 0.030 | 0.005 |
| CYT | +0.018 | −0.002, 0.039 | 0.082 |
| CD8A | +0.013 | −0.005, 0.032 | 0.17 |
| IFNG-6 | +0.012 | −0.009, 0.034 | 0.28 |
| TMB | +0.005 | −0.051, 0.062 | 0.89 |

**Interpretation.**
1. **Compressing the transcriptome into a priori biology is what makes transfer possible.** The same logistic model on all genes drops from 0.609 to 0.555.
2. **The best transferable models sit at about 0.62–0.63 and do not separate from single T-cell markers or TMB.** This matches the field's experience: once the "is there a T-cell infiltrate" signal is captured, bulk RNA adds little that survives across cohorts. I make no claim to have beaten CD8A or IFNG-6.
3. **ImmunoEvo-Net failed, and its ablations show the knowledge mask and adversary did not rescue it.** Its gradient × activation attributions (`program_attribution.csv`) rank KRAS-signalling-down, glycolysis, apical junction, DNA repair and p53 highest. These are **tissue-of-origin and cohort programmes**, not immune biology.
   - The adversary removes cohort identity from the *embedding*. But with 8 training cohorts that each differ in both cancer type and response rate, cohort-correlated biology still predicts labels within the training set.
   - This is a known failure mode of domain-adversarial training when label shift and domain are confounded.
   - The architecture was fixed before evaluation. It was **not** retuned or renamed after these results were seen, and the negative result is reported as is.
4. Per cohort (`loco_auc_by_cohort.csv`, Fig. 2), all methods do well on Gide melanoma (0.80–0.84 for the T-cell scores). Hugo melanoma (n = 26) and Choueiri RCC (n = 16) are close to chance for most methods, and the confidence intervals there are very wide.

## 3. Result 2: tumour-evolution features from clinical panels add little beyond TMB (`evo_*.csv`)

**MSK-IMPACT ICB pan-cancer (n = 1,661; OS; Cox stratified by cancer type, adjusted for drug class, age group, sex and sample type):**
- TMB HR 0.82 per SD.
- With TMB in the model:
  - clonal fraction HR 0.95 (p = 0.31);
  - MATH HR 1.04 (p = 0.35).
- Clonal TMB stays protective (HR 0.74 per SD) when subclonal TMB is added to the model.

5×5 cross-validated C-index:

| Model | C-index |
|---|---|
| Cancer type + drug | 0.617 |
| + TMB | 0.631 |
| + TMB + evolution | 0.633 |
| DeepSurv (+TMB + evolution) | 0.631 |
| RSF (+TMB + evolution) | 0.622 |

**External validation (frozen MSK models).**

| Cohort | Endpoint | C-index, Cox +TMB | C-index, Cox +TMB+evolution |
|---|---|---|---|
| Miao/Van Allen 2018, mixed solid tumours | OS | 0.544 | 0.550 |
| Rizvi 2018, NSCLC | PFS | 0.545 | 0.546 |

- For Rizvi, the 211 patients who also appear in MSK-IMPACT were removed from MSK training before this test.

**Single-feature response AUC:**

| Feature | Allen | Rizvi |
|---|---|---|
| TMB | 0.647 (0.575–0.727) | 0.612 (0.531–0.690) |
| Clonal TMB | 0.637 | 0.578 |
| Clonal fraction | 0.560 | 0.490 |
| MATH (inverse) | 0.560 | 0.470 |

**Why the null?**
- Targeted panels see a few hundred kb, so most tumours have fewer than 10 mutations to estimate clonality from.
- The VAF-based purity estimate agrees only modestly with pathology purity (Spearman 0.40, n = 1,025; `evo_misc.json`).
- The result says that **panel-derived ITH is not a useful ICB biomarker beyond TMB**. It does not say that clonal architecture is irrelevant; exome-based CCF with copy-number correction (as in McGranahan 2016) remains the right test.

## 4. Result 3: what changes in the tumour microenvironment on therapy (`ontx_*.csv`)

**Data:** 77 paired pre/on-treatment biopsies (Riaz 42 pairs, Gide 16, Prince 19; 25 responders).

**Method.**
1. For each programme, compute Δ = on − pre.
2. Compare Δ between responders and non-responders within each cohort with Mann–Whitney.
3. Combine the cohorts with a Stouffer z weighted by √n and apply Benjamini–Hochberg FDR.

**Significant, and same direction in all three cohorts (FDR < 0.05).** These programmes fall more in responders:

| Programme | z |
|---|---|
| mTORC1 | −4.06 |
| Melanocyte pigmentation | −3.64 |
| Glycolysis | −3.59 |
| Unfolded-protein response | −3.39 |
| Hedgehog | −3.00 |
| Protein secretion | −2.93 |
| UV response | −2.88 |
| Mitotic spindle | −2.81 |
| DNA repair | −2.80 |

These are tumour-cell-intrinsic (growth, metabolism and proliferation) programmes. A coordinated drop is what tumour-cell clearance looks like in bulk RNA. CD8 T, exhausted-CD8 and B-cell programmes rise in responders (z ≈ +0.8 to +1.1), but not significantly at this sample size.

## 5. What I would claim, and what I would not

**Claim:**
- In a strict multi-cancer LOCO evaluation, a pre-specified knowledge-programme logistic model with TMB generalises at AUC 0.63. It significantly beats TIDE, GEP-18, PD-L1 and whole-transcriptome learning.
- Deep, knowledge-masked domain-adversarial learning did **not** transfer, and its attributions show why.
- Panel-derived ITH adds nothing beyond TMB.
- On-treatment response coincides with collapse of tumour-intrinsic programmes.

**Not claim:**
- Superiority over CD8A, IFNG-6 or TMB.
- Clinical utility: an AUC of 0.63 is not enough to deny anyone treatment.

## 6. Reproducibility

- **Compute:** all code runs on Google Colab.
- **Data:** cBioPortal datahub tarballs, downloaded by `src/ie.py` on first use.
- **Fixed for all models:** seeds, cohort list and response definition (cBioPortal/iAtlas `RESPONDER` label, pre-treatment samples only).
- **Excluded:** the trained network weights (122 MB) are not in the public repository.
