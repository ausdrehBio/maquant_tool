# Does our workflow reproduce Bergant et al. (2023)? A check of the published proteomics results

*23 November 2026 · Publication: Bergant V, Schnepf D, de Andrade Krätzig N, et al. (2023), Nat Commun 14:4906 · Data: PRIDE PXD040456 · Full tables: `answers.md` in this folder*

## Summary

- **Purpose.** Before using the Bergant et al. proteomics data for new questions, such as the DPR search, we checked whether our own workflow reproduces the published results. The workflow is MaxQuant 2.9.0.0 plus our analysis scripts. We focused on the whole-cell proteome of A549 cells in Fig. 4 and Supplementary Data 5.
- **The fair benchmark.** The published results come from a MaxQuant search that the authors did not deposit. The search they did deposit, with MaxQuant 1.5.6.2, does not match the published table exactly either. We therefore ask whether our workflow comes **as close to the published results as the authors' own deposited search does**.
- **Verdict: 7 of 8 checks consistent, 1 minor deviation, none inconsistent.** On all 21 comparison metrics, our results lie within 5 % of what the authors' deposited search achieves.
- **What we reproduce:**
  - identical search settings
  - all 42 samples assigned to the right condition, with the NS1 mutations visible in the expected strains
  - the same identification depth
  - quantities for the same measurement that differ less than half as much as two replicates do
  - the published statistics, reproduced exactly from the published table
  - fold changes that correlate with the published ones at r = 0.965–0.970
  - the sample structure of Fig. 4a
  - the viral protein pattern of Fig. 4b, and 8 of the 9 tests of the conclusions in Fig. 4c–f.
- **What nobody can reproduce from the deposited data:** the exact number of significant proteins. Both our search and the authors' deposited search find about **30 % fewer** significant proteins than the paper. The published table had far fewer missing values (7 % vs. 23 % filled in by imputation, Section 2). Correspondingly, about half of the individual published hits recur in each re-analysis, but only 2 of 3,592 point in the opposite direction.
- **The minor deviation:** in Fig. 4e, the statistical test reaches p < 0.01 in 80 % of our random imputation draws, against 47 % for the authors' deposited data. The effect itself and its direction are present in every analysis.
- **Conclusion:** the workflow is sound. The DPR search uses the same MaxQuant version and settings, so this reproduction is also the quality check for its basis.

![Overview of the eight checks and their verdicts](figures/summary_verdicts.png)

*Figure 1. The eight checks. Each row shows our value (blue dot) against the verdict rule (green: consistent, yellow: minor deviation, red: inconsistent) and, where the rule uses it, the authors' deposited search (orange diamond).*

---

## 1. What was compared

The study measured the proteome of A549 cells 24 h after infection with six PR8 NS1 mutants (MOI 3) or mock infection, with six replicates per condition: 42 measurements (runs) in total. The NS1 mutants are named after residues 38, 103/106 and 184:

| Condition | NS1 genotype |
|---|---|
| mock | uninfected |
| RSIG | R38, S103/I106, G184: PR8 wild type |
| RSIR | G184R |
| RFMG | S103F/I106M (CPSF30 binding) |
| RFMR | S103F/I106M + G184R |
| AFMG | R38A (dsRNA-binding mutant) + S103F/I106M |
| AFMR | R38A + S103F/I106M + G184R |

The comparisons of "G" with "R" within one background are named RSIX (RSIG vs. RSIR), RFMX (RFMG vs. RFMR) and AFMX (AFMG vs. AFMR).

We compared three versions of the results:

| Version | What it is | MaxQuant | Runs |
|---|---|---|---|
| **Published** | Supplementary Data 5: the final, filtered table with the paper's fold changes and p values. It comes from a search that was **not deposited**. | unknown | 42 |
| **Deposited** | the authors' MaxQuant output, deposited in PRIDE | 1.5.6.2 | 67 (25 not used in the paper) |
| **Ours** | our search of the raw data with the authors' settings | 2.9.0.0 | the paper's 42 |

Two pieces of evidence show that the published table comes from a different search than the deposited one. First, it uses protein identifiers from a newer Ensembl release: 72 of its protein groups do not exist in the deposited sequence database. Second, it has far fewer missing values. No re-analysis of the deposited data can therefore match the published table exactly. The authors' deposited search shows how close one can realistically get, and we use it as the yardstick throughout.

## 2. Background: what you need to know about quantitative proteomics

How proteomics identifies peptides, and how false identifications are controlled, is explained in Sections 2–3 of the DPR report (`../../Bergant_2023_adapted/analysis_of_results/DPR_evidence_report.md`). For this report, five further concepts matter:

- **Protein groups.** MaxQuant identifies peptides and assembles them into proteins. Peptides shared by several proteins (isoforms, paralogues) cannot be assigned unambiguously, so MaxQuant reports a **protein group**: the set of proteins that the peptides found cannot tell apart. Two searches can draw these groups slightly differently. We therefore match groups between tables through shared protein identifiers, and count groups that were split or merged separately.
- **Label-free quantification (LFQ).** Without labels, the amount of a protein in a sample is estimated from the signal strength (**intensity**) of its peptides in the mass spectrometer. MaxQuant's **LFQ** algorithm normalises these intensities so that the same protein can be compared across samples. All values are on a log2 scale, so a difference of 1 means a two-fold change.
- **Match between runs.** If a peptide was sequenced in one run but not in another, MaxQuant can still find its signal in the second run from its mass and retention time. This fills gaps but is less certain than a direct identification.
- **Missing values and imputation.** A protein near the detection limit is often not detected in every sample. The paper, like most label-free studies, first removes proteins that are not detected in at least 5 of 6 replicates of at least one condition. It then fills the remaining gaps with random values drawn from a distribution at the low end of all measured values. This is called **imputation**, and its parameters are given in the Methods. Because imputation is random, results near the threshold can change from one draw to the next. The authors published one draw. We repeated the imputation **100 times** with different random seeds and report medians and ranges.
- **Statistics.** For each protein, the paper fits a linear model and tests each strain against mock (t test) and each G strain against its R counterpart (z test). The p values are corrected for multiple testing (Benjamini–Hochberg). A protein counts as significant if it changes more than 1.7-fold (|log2 fold change| > 0.75) with an adjusted p < 0.01. Note that "false discovery rate" is used in two senses: the decoy-based error rate of identification (DPR report) and the multiple-testing correction of p values used here.

## 3. The eight checks

| | Question | Published / expected | Our result | Verdict |
|---|---|---|---|---|
| Q1 | Do our search settings match the authors'? | MaxQuant 1.5.6.2, default parameters, LFQ, match between runs | 0 relevant differences in 39 settings | consistent |
| Q2 | Are the samples assigned to the right conditions? | 42 runs, 6 per condition | 42/42 correct; NS1 mutations visible as expected | consistent |
| Q3 | Do we identify as much? | 6,052 protein groups after the paper's filter | 1.03× the authors' identifications per run; 95.8 % of published groups found | consistent |
| Q4 | Are the quantities the same for the same measurement? | (deposited values) | difference 0.44× that between two replicates | consistent |
| Q5 | Does our code reproduce the published statistics? | Supplementary Data 5 | exactly, every call identical | consistent |
| Q6 | Do the differential results agree? | e.g. RSIG vs. mock: 124 up / 205 down | fold changes r = 0.965–0.970; counts as in the authors' deposited search | consistent |
| Q7 | Is the sample structure the same? | Fig. 4a: CPSF30-binding strains separate from G184R strains | reproduced (similarity 0.986) | consistent |
| Q8 | Do the paper's conclusions hold? | Fig. 4b–f | viral pattern reproduced (r = 0.91); 8 of 9 claim tests hold | minor deviation |

### Q1. Search settings: identical

We compared 39 settings that influence identification or quantification. These include:
- enzyme, missed cleavages and modifications
- mass tolerances and error rates
- match between runs and the LFQ parameters.

**None differs.** The only intended differences are the MaxQuant version and the set of runs (the paper's 42 instead of all 67). Five further settings differ, all belonging to features this experiment does not use (de novo sequencing, cross-linking). Every protein identifier in our results exists in the authors' sequence database (human Ensembl GRCh37, PR8, SC35M).

### Q2. Sample assignment: all 42 samples correct

Raw data files have cryptic names, and several samples were measured a second time. We therefore checked independently that each file belongs to the condition it is assigned to.

- **Protein profiles.** Every published sample resembles our sample of the same name more than any other sample (42/42). The margin to the most similar sample of another condition is comfortable (difference in correlation ≥ 0.40).
- **NS1 genotype.** The NS1 mutations can be read directly from the data:
  - The wild-type G184 peptide (`NAVGVLIGGLEWNDNTVR`) was sequenced in all six samples of each G strain (RSIG, RFMG, AFMG) and in none of the R strains or mock.
  - The wild-type S103/I106 peptide (`DWSMLIPK`) was found only in RSIG and RSIR.
  - Residue 38 cannot be tested, because its tryptic peptides are too short.
- **Viral traces in mock.** Mock samples contain a tiny viral signal (3 viral spectra vs. 32,895 in infected samples, 0.01 %). The authors' deposited search shows the same. It most likely reflects carryover between measurements on the instrument, not a sample mix-up.

### Q3. Identifications: the same depth

| Per run (median of 42) | Authors' deposited search | Ours | Ratio |
|---|---|---|---|
| protein groups identified by MS/MS | 4,693 | 4,850 | 1.03 |
| MS/MS spectra assigned | 39,996 | 40,778 | 1.02 |
| protein groups added by match between runs | 1,122 | 791 | 0.70 |
| protein groups with an LFQ value | 4,631 | 4,440 | 0.96 |

Our search identifies slightly more by direct sequencing but transfers fewer proteins by match between runs. Two plausible causes exist: the stricter matching algorithm of the newer MaxQuant, and the smaller number of runs available to transfer from (42 instead of 67). A search of all 67 runs with the new version would separate the two; it is set up but not yet analysed.

Of the 6,052 protein groups in the published table, our search finds 95.8 % and the authors' deposited search 97.3 %. The published groups we miss are mostly low-abundance proteins near the detection limit. The other search found them and ours did not, which is expected near the limit.

### Q4. Quantities: closer than replicates

For the same raw file, we compared our LFQ value of each protein with the authors' deposited value. Our values differ from theirs by a typical 0.16 log2 units. Two replicate samples of the same condition differ by 0.35. **The difference between the two software versions is less than half the natural variation between replicates.** Neither table is closer to the published values than the other (deviation 0.211 for ours, 0.210 for the deposited search).

Our LFQ values are lower overall by a constant amount. This is how LFQ scales intensities within each search, and a constant shift does not change any fold change.

### Q5. Statistics: reproduced exactly

Starting from the published table, our code reproduces **every published fold change, p value and significance call exactly**: the largest deviation is about 10⁻¹². The Methods describe the statistics only briefly. Reproducing them exactly established which tests the authors used: a linear-model t test for strain vs. mock, a z test for G vs. R, and Benjamini–Hochberg correction within each comparison.

We could also tell which published values are real measurements and which were imputed. MaxQuant writes intensities with five significant digits, and imputed values lack this rounding. In the published table, 92.6 % of the values are real. The imputed ones have exactly the mean and spread that the published imputation parameters predict. When we re-impute the published table 100 times with our code, 16 of the 18 published up/down counts fall within our range of outcomes. The other two (AFMG vs. mock, AFMR vs. mock) lie just outside, a normal consequence of the authors' single random draw.

### Q6. Differential results: as close as the authors' own data

We ran the paper's pipeline on our search and on the authors' deposited search, and compared both with the published results.

- **Fold changes.** For proteins without imputed values, our fold changes correlate with the published ones at **r = 0.965–0.970** for the six strain vs. mock comparisons. The deposited search gives 0.966–0.974. For the G vs. R comparisons the values are 0.915–0.963 vs. 0.932–0.960.
- **Direction.** Among the published significant proteins, 92–96 % change in the same direction in our data.
- **Counts.** Both re-analyses find about 30 % fewer significant proteins than the paper, and to the same extent:

| Comparison | Published (up / down) | Authors' deposited search | Ours (median of 100 draws) |
|---|---|---|---|
| RSIG vs. mock | 124 / 205 | 83 / 149 | 80 / 153 |
| RSIR vs. mock | 174 / 247 | 112 / 186 | 113 / 178 |
| RFMG vs. mock | 85 / 371 | 59 / 260 | 54 / 255 |
| RFMR vs. mock | 184 / 318 | 124 / 217 | 121 / 223 |
| AFMG vs. mock | 110 / 323 | 65 / 225 | 61 / 238 |
| AFMR vs. mock | 229 / 246 | 175 / 156 | 170 / 149 |
| RSIX (RSIG vs. RSIR) | 27 / 65 | 14 / 47 | 13 / 49 |
| RFMX (RFMG vs. RFMR) | 45 / 294 | 41 / 204 | 44 / 182 |
| AFMX (AFMG vs. AFMR) | 141 / 404 | 95 / 300 | 85 / 315 |

![Significant proteins per comparison: published, deposited and ours](figures/q06_counts.png)

*Figure 2. Significant proteins per comparison (up above the line, down below) in the published table (black), the authors' deposited search (orange) and our search (blue).*

The reason is the number of missing values. Before imputation, the published table had real values for 92.6 % of its entries, the deposited search for 76.9 % and ours for 77.0 %. Imputed values are random and dilute real differences, so fewer proteins pass the threshold. Whatever produced the published table delivered more complete data than either MaxQuant search of the deposited files.

On average we call 54 % of the individual published hits; the deposited search calls 55 %. Most of the published hits that are lost involve imputed values or sit just below the threshold. **Only 2 of the 3,592 published calls change direction in our data.** Disagreement is concentrated in low-abundance proteins and in proteins with many imputed values, as expected.

### Q7. Sample structure (Fig. 4a): reproduced

Fig. 4a shows a map of the samples (UMAP) in which the CPSF30-binding strains lie apart from the G184R strains. UMAP maps are random and change from run to run. We therefore compared the sample structure with principal component analysis (PCA), a deterministic method that summarises all protein values into a few axes. We then measured how well our sample map can be superimposed on the published one, using Procrustes similarity, where 1 means identical.

Our map matches the published one with a similarity of **0.986** (deposited search: 0.991). The separation shown in Fig. 4a holds in our data: the distance between the CPSF30-binding and the G184R strains is larger than the distances within either group (40.1 vs. 30.8; published 48.2 vs. 35.9).

### Q8. The conclusions of Fig. 4: 8 of 9 tests hold

We translated the statements of Fig. 4b–f into fixed tests and ran them on each of our 100 imputation draws. A conclusion counts as reproduced if it holds in at least 95 % of the draws.

| Fig. | Statement | Published table, re-imputed | Authors' deposited search | Ours | Verdict |
|---|---|---|---|---|---|
| 4b | pattern of viral protein levels across the strains | – | r = 0.945 | r = 0.910 | consistent |
| 4c | fewer G vs. R differences in the wild-type background (RSIX) than in the CPSF30-binding backgrounds (RFMX, AFMX) | 100 % | 100 % | 100 % | consistent |
| 4e | in the CPSF30-binding backgrounds, fold changes of G strains are shifted downward relative to R strains (statistical test, p < 0.01) | 100 % | 47 % | 80 % | **minor deviation** |
| 4e | the share of strongly decreased proteins is larger in RFMX and AFMX than in RSIX | 100 % | 100 % | 100 % | consistent |
| 4f | antiviral proteins (GO: defence response to virus) are induced by the wild-type virus | 100 % | 99 % | 99 % | consistent |
| 4f | RSIR resembles RSIG for these proteins | 100 % | 100 % | 99 % | consistent |
| 4f | they are lower in RFMG than in RSIG | 100 % | 100 % | 100 % | consistent |
| 4f | they are lower in AFMG than in RSIG | 100 % | 100 % | 100 % | consistent |
| 4f | G184R relieves the suppression in both CPSF30-binding backgrounds (RFMX, AFMX) | 100 % | 100 % | 100 % | consistent |

Percentages give the share of the 100 imputation draws in which the statement holds. Apart from 4b, which is a correlation, the table contains nine tests (the last row combines two, one per background), and eight of them hold.

![Fold-change distributions and antiviral proteins, published vs ours](figures/q08_claims.png)

*Figure 3. Top row: published table; bottom row: ours. The first three panels show the cumulative distributions of fold changes vs. mock for the G and R strain of each background (Fig. 4e). In the RFM and AFM backgrounds the G strain (dark) lies to the left of the R strain (light), meaning more decreased proteins. The last panel shows the fold changes of the antiviral GO set (Fig. 4f).*

**The minor deviation (Fig. 4e).** The downward shift of the G strains is present in our data, with the published direction. In the median of our draws the test reaches p = 0.002 in the RFM and p = 0.005 in the AFM background. Only in 80 % of the draws, however, does it stay below 0.01 in *both* backgrounds at once. The authors' own deposited data fare worse (47 %). The reason is the same as in Q6: more imputed values (23 % instead of 7 %) blur the distributions. One further statement, that the shift is stronger in the CPSF30-binding backgrounds than in the wild-type background, holds in only 86 % of re-imputations of the *published* table itself. It therefore depends on the authors' single random draw and was reported but not judged.

## 4. Observations along the way

- **NS1 of the AFMR virus.** The wild-type NS1 peptide covering residues 22–37 (`VADQELGDAPFLDRLR`) was sequenced in all six samples of five strains but in **none of the AFMR samples**. Other NS1 peptides are found in AFMR, and AFMG, which carries the same R38A mutation, has this peptide. This may point to an additional NS1 sequence difference near residues 22–38 in the AFMR virus that its name does not describe. The authors' NS1 sequences would be needed to tell. The DPR analysis also flagged an unusual peptide from this region (`VADQELGDAPFLDRR`).
- **NS1 quantities.** The sequence database contains only wild-type PR8 NS1. Peptides covering mutated residues are therefore missed in some strains, so differences in NS1 levels between strains partly reflect which peptides can be detected, not only how much NS1 is present.
- **PB2 grouping.** The new MaxQuant version reports PR8 PB2 together with its SC35M counterpart in one protein group, whereas the old version kept them apart. Both are in the database because the authors searched against both viruses.

## 5. Limitations

- The published table comes from a search that is not deposited, so the published significance counts cannot be matched by any re-analysis of the deposited data. We used the authors' deposited search as the realistic benchmark.
- Our search used the paper's 42 runs, the deposited search all 67. Until the 67-run search with the new version has been analysed, the effect of the software version cannot be separated from the effect of the run set. This matters mainly for match between runs (Q3).
- Imputation is random. We report medians and ranges over 100 draws. Published conclusions that depend on the authors' single draw are reported but not judged.
- The antiviral protein set (GO:0051607) and the gene names come from a recent annotation snapshot, not from the authors' annotation version. This explains the slightly different set sizes (89 vs. 94 protein groups).

## 6. What this means for the DPR analysis

The DPR search (`Bergant_2023_adapted`) uses the same MaxQuant version and the same settings as the search checked here. It differs in two ways: it includes the 25 additional deposited runs (A149V, SC35M, SC35MdNS1, a second mock series and one repeat measurement), and it adds the DPR sequences to the database. This reproduction shows that the underlying workflow is sound:
- the settings match the authors'
- the samples are correctly assigned
- identification depth and quantification are equivalent to the authors' own processing
- the paper's biological conclusions hold.

The problems described in the DPR report, such as the distortion of the viral protein groups, are caused by the added DPR sequences, not by the workflow.

---

## Glossary

| Term | Meaning |
|---|---|
| Benjamini–Hochberg correction | adjustment of p values for the thousands of proteins tested at once; the "adjusted p" |
| Deposited search | the authors' own MaxQuant output, published with the raw data in PRIDE |
| Fold change (log2) | difference between two conditions on a log2 scale; 1 = two-fold, 0.75 ≈ 1.7-fold |
| Imputation | filling in missing values with random low values before statistical testing |
| Intensity | the signal strength of a protein's peptides, a measure of its amount |
| LFQ | MaxQuant's label-free quantification: intensities normalised for comparison across samples |
| Match between runs | finding a peptide's signal in a run where it was not sequenced, from its mass and retention time |
| PCA | principal component analysis: summarises thousands of protein values into a few axes to map samples |
| Procrustes similarity | how well two sample maps can be superimposed after shifting, rotating and scaling (1 = identical) |
| Protein group | proteins that cannot be told apart with the peptides found, reported together |
| Replicate | an independent sample of the same condition (here six per condition) |
| Run | one mass-spectrometry measurement of one sample |
| Seed (imputation draw) | one repetition of the random imputation; we used 100 |
| UMAP | a non-linear method for mapping samples, used in Fig. 4a; its result changes with each run |

## Files

- `answers.md` contains the complete, automatically generated results of the eight checks, including each verdict rule and all intermediate numbers.
- `figures/` and `tables/` contain all figures and per-protein results (for example, our fold changes and significance calls for every protein).
- `scripts/` contains the analysis. `mq_common.py` documents all definitions, and `run_all.py` runs the checks.
