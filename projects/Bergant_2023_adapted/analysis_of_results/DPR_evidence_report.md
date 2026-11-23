# Are DPRs detectable in the proteome? A re-analysis of the Bergant et al. (2023) influenza proteomics data

*22 November 2026 · Data: PRIDE PXD040456 · Search: MaxQuant 2.9.0.0 · Full tables: `answers.md` in this folder*

## Summary

- We re-searched the published mass-spectrometry data of A549 cells (67 measurements: mock, seven PR8 variants, H7N7 SC35M) against the usual human and influenza proteins **plus 58,397 DPR sequences**. We then asked which peptides can only have come from a DPR.
- At MaxQuant's standard error rate of 1 %, the search reports **317 such peptides**. For this particular group, however, the share of false hits is about **seven times higher** than for the search as a whole (2.5 % vs. 0.4 %). The standard result is therefore not reliable on its own.
- Nine additional filters, explained below, reduce the list to **37 peptides**. About **one** of them is expected to be a random false hit; statistically, it could be up to about five.
- The core evidence is **24 peptides**: **8** span a deletion junction and **16** run from viral sequence into a new reading frame. Fifteen come from NP, four from M1, three from HA and two from PB1. Of the other 13 peptides:
  - eight differ from a viral protein by a single residue that sits at a junction
  - four contain no recognisable viral sequence
  - one turned out on closer inspection to be the ordinary PB1-F2 of the virus strain.
- Nearly every DPR is supported by **a single peptide** (there are two exceptions). **No spectrum has yet been inspected or validated**, and three of the 37 peptides also appear in H7N7-infected cells. The 37 are therefore **candidates**, not confirmed detections.
- There is a side effect. Adding the DPR sequences distorts how the ordinary viral proteins are reported: M1 disappears and NS1 loses about 98 % of its spectra. This search should be used for the DPR question only.
- Next steps (Section 9):
  - inspect the spectra of the core peptides
  - check that their junctions match those found by DelVG sequencing
  - confirm the best candidates with synthetic peptides.

---

## 1. Data and approach

The data are the whole-cell proteomics measurements of Bergant et al. (2023), *Nat Commun* 14:4906, deposited in PRIDE as PXD040456. Each sample of A549 cells was measured by mass spectrometry, and one measurement is called a **run**. The paper used 42 runs (mock and the six NS1 variants). The other 25 runs (A149V, SC35M, SC35MdNS1, a second mock series and one repeat measurement) are deposited with the data but not used in the paper. We included all 67.

| Condition | Virus | Runs |
|---|---|---|
| mock, mock2 | uninfected (two series) | 6 + 6 |
| RSIG | PR8 wild type (NS1 R38, S103/I106, G184) | 7 (one sample measured twice) |
| RSIR | PR8, NS1 G184R | 6 |
| RFMG | PR8, NS1 S103F/I106M (CPSF30 binding) | 6 |
| RFMR | PR8, NS1 S103F/I106M + G184R | 6 |
| AFMG | PR8, NS1 R38A + S103F/I106M | 6 |
| AFMR | PR8, NS1 R38A + S103F/I106M + G184R | 6 |
| A149V | PR8, NS1 A149V | 6 |
| SC35M | H7N7 SC35M | 6 |
| SC35MdNS1 | H7N7 SC35M lacking NS1 | 6 |

We processed the data with **MaxQuant**, the standard software for this type of measurement. We used the same settings as in our reproduction of the published analysis, with one change: the DPR sequences (`PR8_3frame_uniq_*`, 58,397 entries) were added to the reference sequences. The reference sequences are the human proteome (Ensembl GRCh37, 104,763 entries), the PR8 and SC35M proteins, and a standard list of laboratory contaminants such as keratins and trypsin. An analysis script then extracted every peptide that occurs only in DPR sequences and filtered these peptides step by step (Section 5).

## 2. How proteomics identifies peptides

### 2.1 From protein to spectrum

The method used here, "bottom-up" proteomics, does not measure whole proteins. It measures short pieces of them.

1. **Digestion.** The proteins are cut with the enzyme trypsin, which cleaves after lysine (K) and arginine (R). The resulting **peptides** are typically 7–30 residues long and, with the exception of protein C-termini, all end in K or R. Trypsin occasionally skips a site; the search allowed up to two such *missed cleavages* per peptide.
2. **Separation.** The peptide mixture flows through a chromatography column, and different peptides leave it at different times. The time at which a peptide leaves is its **retention time**, which is reproducible for a given peptide.
3. **Mass measurement (MS1).** The peptides leaving the column are given an electric charge, usually two to four protons. The mass spectrometer measures their **m/z**, the mass divided by the charge. For example, a peptide of 1,500 Da carrying two charges appears at m/z ≈ 751. The measurement is extremely precise. MaxQuant accepted deviations of at most 4.5 **ppm** (parts per million), which for a 1,500 Da peptide is 0.007 Da, roughly 1/140 of the mass of a hydrogen atom.
4. **Fragmentation (MS/MS).** The instrument isolates one peptide ion at a time, breaks it apart and measures the masses of the fragments. The result is an **MS/MS spectrum**. Peptides break mainly at the peptide bonds, which produces two "ladders" of fragments, one from each end. The mass difference between neighbouring rungs equals the mass of one amino acid. An MS/MS spectrum is therefore a fingerprint of the sequence. In practice the fingerprint is usually incomplete, because some rungs are missing or weak.

Each run produces tens of thousands of MS/MS spectra.

### 2.2 From spectrum to sequence: the database search

MaxQuant does not read sequences directly from the spectra. Instead, it compares each spectrum with a list of candidate sequences, the **database**: the protein sequences we supply. It digests every protein in the database *in silico*. For each spectrum it takes the candidate peptides whose mass fits, predicts their fragment ladders and scores how well each prediction matches the measured spectrum. The best-scoring candidate becomes the identification, a **peptide-spectrum match**.

Three consequences matter for the DPR question:

- **Only sequences that are in the database can be found.** This is why the DPR sequences had to be added. Without them, DPR peptides would be invisible.
- **Every spectrum is assigned to its best-fitting candidate, even when the true peptide is not in the database.** Some spectra come from peptides missing from the database, for example from a human sequence variant or a chemically modified peptide. Such a spectrum can be assigned to a DPR sequence that happens to fit reasonably well. This is the main source of false DPR identifications.
- **A larger database produces more chance matches.** The DPR sequences enlarge the database by more than half (58,397 entries on top of 104,763 human proteins).

### 2.3 What mass spectrometry cannot tell apart

- **Isoleucine and leucine** have identical masses and are treated as the same residue.
- **Combinations of residues** can match the mass of a single residue. Glycine + valine, for example, weigh the same as arginine. A peptide containing "GV" can therefore masquerade as one containing "R" (see filter 6).
- **Chemical modifications.** Sample preparation changes some peptides chemically, and the search has to allow for this. Here MaxQuant allowed three modifications:
  - oxidation of methionine, an artefact of sample handling
  - acetylation of protein N-termini, which is biological
  - **deamidation** of asparagine (N) and glutamine (Q), which happens during sample preparation.

  Deamidation turns N into aspartate (D) and Q into glutamate (E). The product is chemically identical to an encoded D or E. A peptide that differs from a known protein only by N→D or Q→E is therefore indistinguishable from the known peptide after deamidation.

## 3. How false identifications are controlled

### 3.1 Decoys: a negative control built into the search

MaxQuant adds a reversed copy of every protein to the database, called a **decoy**. Decoys have the same composition and length as real proteins but cannot occur in the sample. Every match to a decoy is therefore false by construction. False matches to real sequences (the **targets**) arise by chance at about the same rate. Counting decoy hits therefore estimates how many target hits are false.

### 3.2 False discovery rate (FDR) and posterior error probability (PEP)

- **FDR.** MaxQuant ranks all identifications by confidence and accepts them down to the point where the decoys make up 1 % of the accepted list. This is the standard "1 % FDR": about one in a hundred accepted identifications is expected to be false. The FDR describes the list as a whole, not any single identification in it.
- **PEP.** For each identification, the PEP is MaxQuant's estimate of the probability that this particular identification is wrong. A list at 1 % FDR can contain individual identifications with a PEP of several percent, as long as the list as a whole stays at 1 %. In our DPR list, 74 of the 317 peptides have a best PEP above 5 %, and the highest is almost 9 %. "PEP < 0.01" means a less than 1 % chance of being wrong.

### 3.3 Why "1 % overall" does not mean "1 % for DPR peptides"

The FDR is controlled for the entire list, but false hits are not evenly distributed across it. They make up a larger share wherever true hits are rare. A diagnostic test with a fixed false-positive rate behaves the same way: applied to a rare condition, most of its positive results are false.

The same thing happens here. Most of the roughly 90,000 peptides identified are ordinary human peptides, a group with abundant true hits. DPR-specific peptides are few and new, and they are drawn from a large added search space, so false hits make up a much larger share of them:

| | Peptides | Decoys | Decoys / peptides |
|---|---|---|---|
| Whole search | ~89,800 | 343 | 0.4 % |
| DPR-specific peptides only | 317 | 8 | 2.5 % |

Decoys can be classified as "DPR-specific" in exactly the same way as real peptides. We can therefore follow the decoys through every filter, and the decoys that remain estimate the false hits that remain (Section 5).

## 4. Which peptides count as DPR-specific

A peptide counts as DPR-specific if all three of the following hold:

1. MaxQuant identified it at 1 % FDR.
2. It is a tryptic peptide of at least one DPR sequence.
3. It occurs nowhere in the human, PR8, SC35M, contaminant or control sequences.

The same definition applied to the reversed sequences gives the DPR-specific decoys.

To make the peptides interpretable, each one was compared with the PR8 proteins and placed in one of four classes. The classes describe the peptides; they are not filters.

| Class | Definition | Interpretation |
|---|---|---|
| **Junction** | begins in one region of a viral protein and continues, in the same reading frame, with a distant region (each side ≥ 5 residues) | spans an in-frame deletion junction; the most direct evidence |
| **Frame continuation** | one end matches ≥ 5 residues of a viral protein; the rest is found in no viral protein (the new part may precede or follow the viral part) | crosses a junction into another reading frame, or extends beyond the normal end of the protein |
| **Substitution** | differs from a viral protein in exactly one residue | a junction with only one residue beyond it, a sequence variant of the virus, or a modification (Section 6.3) |
| **No viral anchor** | neither end matches ≥ 5 residues of a viral protein | an entirely out-of-frame translation, or a human peptide missing from the database |

**Notation in the tables:**
- "NP 33-39 → NP 420-436": the first part of the peptide equals NP residues 33–39 and the second part equals NP residues 420–436.
- "3 new + NP 279-293": three residues found in no viral protein, followed by NP residues 279–293.
- "NS1 L36R": identical to NS1 except that L36 is replaced by R.

## 5. Filters and what each one removes

The filters are applied one after another, and each step keeps only the peptides that pass every filter so far. Decoys go through exactly the same steps.

| # | Filter | Requirement | Rationale |
|---|---|---|---|
| 1 | MS/MS | sequenced by at least one MS/MS spectrum | MaxQuant can also transfer an identification from one run to another based on mass and retention time alone (*match between runs*), without any fragment spectrum in that run. A peptide known only from such transfers was never actually sequenced. |
| 2 | Main search | at least one identification from the regular search | MaxQuant also searches spectra that already have an identification for a second, weaker peptide hidden in the same spectrum (*second-peptide search*). These identifications are more error-prone. |
| 3 | PEP < 0.01 | best identification has an error probability below 1 % | The 1 % FDR admits individual identifications that are far less certain (Section 3.2). |
| 4 | Not in mock | never sequenced in uninfected cells (mock, mock2) | A virus-derived peptide cannot be genuine in uninfected cells. |
| 5 | Not SC35M-only | not sequenced exclusively in H7N7-infected cells | The DPR sequences derive from PR8. |
| 6 | No isobaric twin | no known sequence of the same mass (within 10 ppm) that shares at least half of the peptide | The spectrum could equally come from that known peptide. Example: `RQIASNENMETMESSTLELR` (815 spectra) has the same mass as the ordinary NP peptide `GVQLASNENMETMESSTLELR`, because G + V weigh the same as R. |
| 7 | Not a deamidated known peptide | at least one identification whose deamidated residues, written as D or E, do not reproduce a known sequence | Otherwise the "DPR" spectrum is simply a known peptide that was deamidated during sample preparation (Section 2.3). |
| 8 | ≥ 3 runs | regular identifications in at least 3 of the 67 runs | Reproducibility: a genuine peptide should be found in more than one or two samples. |
| 9 | One spectrum without deamidation | at least one identification that needs no deamidation to fit | The strictest step: the new residues must fit the spectrum without assuming any chemical change. |

| Step | After filter | Peptides left | Removed | Decoys left | Decoys / peptides | Junction | Frame cont. | Substitution | No anchor |
|---|---|---|---|---|---|---|---|---|---|
| L0 | MaxQuant 1 % FDR | 317 | – | 8 | 2.5 % | 64 | 166 | 51 | 36 |
| L1 | MS/MS | 314 | 3 | 8 | 2.5 % | 64 | 163 | 51 | 36 |
| L2 | Main search | 273 | 41 | 8 | 2.9 % | 51 | 144 | 45 | 33 |
| L3 | PEP < 0.01 | 110 | 163 | 7 | 6.4 % | 24 | 54 | 22 | 10 |
| L4 | Not in mock | 104 | 6 | 4 | 3.8 % | 24 | 52 | 22 | 6 |
| L5 | Not SC35M-only | 99 | 5 | 4 | 4.0 % | 24 | 48 | 21 | 6 |
| L6 | No isobaric twin | 96 | 3 | 4 | 4.2 % | 24 | 46 | 20 | 6 |
| L7 | Not a deamidated known peptide | 95 | 1 | 4 | 4.2 % | 24 | 46 | 19 | 6 |
| L8 | ≥ 3 runs | 51 | 44 | 1 | 2.0 % | 12 | 23 | 12 | 4 |
| L9 | One spectrum without deamidation | **37** | 14 | **1** | 2.7 % | **8** | **16** | **9** | **4** |

![Figure 1: number of DPR-specific peptides left after each filter, by class](figures/dpr_filter_cascade.png)

*Figure 1. DPR-specific peptides left after each filter, coloured by class. The numbers in brackets give the decoys left and the DPR sequences that contain at least one of the peptides.*

Three features of the cascade are worth pointing out:

- **The PEP filter removes the most peptides (163), but not the most decoys.** Only one of eight decoys fails it, so the decoy share actually *rises*, from 2.9 % to 6.4 %. The decoys that survive have good scores because their spectra come from real, abundant peptides that are missing from the database, and the search assigned those spectra to reversed sequences. Three of these decoys were identified reproducibly in 35–49 runs, including mock runs. The mock filter removes them; the PEP filter does not. Many confident spectra across many runs are therefore no proof on their own.
- **The mock and reproducibility filters remove six of the eight decoys**, three each.
- **At the end, 1 decoy remains among 37 peptides (2.7 %).** With 37 peptides, 2.7 % is the smallest non-zero value possible. The estimate cannot be more precise than "about one false hit".

## 6. Results

### 6.1 Overview

Thirty-seven peptides pass all filters: 8 junction, 16 frame continuation, 9 substitution and 4 without a viral anchor. The 24 junction and frame-continuation peptides form the core evidence. Fifteen of them derive from NP, four from M1, three from HA and two from PB1.

Proteomics detects abundant proteins far more readily than rare ones, and some sequences yield peptides that are intrinsically hard to detect. That no core peptide was found for a given segment (PB2, PA, NA, NS) is therefore not evidence that its DPRs are absent.

Several downstream regions recur, each joined to different upstream sequences:
- NP 244–261: `…RNPGNAEFEDLTFLAR`, three peptides
- NP 369–382: `…NMETMESSTLELR`, three peptides
- NP 420–436: `…MAAFNGNTEGR`, three peptides
- M1 187–210 and 227–242: the same non-M1 stretch `…GTGIQIT…` joined to different M1 positions.

This could reflect breakpoint hotspots in the DelVG population. It also means that the shared, ordinary viral part may dominate these spectra, so the evidence for the junction rests on the few fragments that cover it (Section 7).

### 6.2 The 37 peptides

**Columns:**
- **Best PEP:** error probability of the best identification. Smaller is better; 1.1e-15 means 1.1 × 10⁻¹⁵.
- **Spectra:** number of regular identifications across all runs.
- **Runs:** number of runs with at least one such identification.
- **Conditions:** conditions in which the peptide was sequenced at least once.

Within each class, peptides are sorted by the number of runs.

#### Junction (8)

| Peptide | Viral origin | Best PEP | Spectra | Runs | Conditions | Note |
|---|---|---|---|---|---|---|
| `EEIRRENMETMESSTLELR` | NP 114-118 (or PB1 751-755) → NP 369-382 | 7.5e-09 | 11 | 7 | RFMR, SC35M | also in 2 SC35M runs |
| `QTYDWTLNDNMTK` | PB1 127-134 → PB1 193-197 | 5.8e-03 | 7 | 7 | RSIG, RSIR, RFMR, AFMR, SC35M | also in 3 SC35M runs |
| `MCTELKLSDYEGRLIENMETMESSTLELR` | NP 43-57 → NP 369-382 | 2.8e-04 | 5 | 5 | RSIG, AFMG, A149V | |
| `CNIAGFLDIWTYNAELLVLLENER` | HA 72-76 → HA 430-449 | 1.1e-15 | 4 | 4 | RSIG, RSIR, RFMG, RFMR, AFMG, AFMR, A149V | |
| `IGGIGRFDRTTIMAAFNGNTEGR` | NP 33-39 → NP 420-436 | 3.1e-03 | 4 | 4 | RSIR, AFMG, A149V | |
| `MVTTTNPLIKAMEQMAGSSEQAAEAMEVASQAR` | M1 165-173 → M1 187-210 | 1.2e-08 | 3 | 3 | RSIR, RFMR | |
| `VLSAFDERRGQISIQPTFSVQR` | NP 67-75 → NP 404-416 | 1.2e-06 | 3 | 3 | AFMG, A149V | |
| `IGGIGRFYIQMCTELMAAFNGNTEGR` | NP 33-47 → NP 425-436 | 7.5e-04 | 3 | 3 | RSIG, RSIR, AFMG | |

#### Frame continuation (16)

| Peptide | Viral origin | Best PEP | Spectra | Runs | Conditions | Note |
|---|---|---|---|---|---|---|
| `MQLCVYGPAVASGYDFER` | 3 new + NP 279-293 | 7.6e-05 | 41 | 33 | RSIG, RSIR, RFMG, RFMR, AFMG, AFMR, A149V | |
| `RMENFLDIWTYNAELLVLLENER` | 4 new + HA 431-449 | 1.1e-21 | 20 | 20 | RSIG, RSIR, RFMG, RFMR, AFMG, A149V | |
| `LIQNSLTIESDKEEIR` | NP 56-64 + 7 new | 3.1e-03 | 17 | 16 | RSIG, RSIR, RFMG, RFMR, AFMG, AFMR, A149V | |
| `GTGIQITWTKQLNCKAMEQMAGSSEQAAEAMEVASQAR` | 14 new + M1 187-210 | 1.1e-08 | 18 | 14 | RSIG, RSIR, RFMG, RFMR, AFMG, A149V | |
| `RSGGVGTMVMELVR` | 3 new + NP 185-195 | 1.2e-13 | 13 | 13 | RSIG, RSIR, RFMG, AFMG, A149V | |
| `FYIQMCLRSR` | NP 39-44 + 4 new | 2.2e-03 | 10 | 10 | RSIR, RFMG, RFMR, AFMR, SC35M | also in 1 SC35M run |
| `LLQNSQVYSLIRPNDGIKYT` | NP 306-319 + 6 new | 3.4e-04 | 10 | 7 | RSIG, RSIR, RFMG, RFMR, AFMG, A149V | |
| `MNRDQAINGITNK` | 5 new + HA 387-394 | 9.7e-04 | 7 | 7 | RSIG, RSIR, RFMG, AFMG | |
| `GTGIQITWTKQLNCAGLKNDLLENLQAYQK` | 14 new + M1 227-242 | 5.0e-03 | 8 | 7 | RSIG, RSIR, RFMG, RFMR, AFMR | |
| `GEKESRNPGNAEFEDLTFLAR` | 3 new + NP 244-261 | 1.4e-03 | 5 | 5 | RSIG, RSIR, RFMG, A149V | |
| `QLITGMDPRMCSLMQGSTLPR` | 3 new + NP 157-174 | 2.2e-03 | 6 | 5 | RFMG, RFMR, AFMG, A149V | |
| `GRWTTNTETEK` | PB1 53-61 + 2 new | 2.3e-06 | 5 | 4 | A149V | |
| `MPLMGTGIQITKNDLLENLQAYQK` | 11 new + M1 230-242 | 4.7e-04 | 4 | 4 | RSIR, RFMG, RFMR, A149V | |
| `DEKQELLMSRNPGNAEFEDLTFLAR` | 8 new + NP 245-261 | 6.9e-04 | 4 | 4 | RSIR, RFMR, AFMG, A149V | |
| `NKMMESARPEDVSFQGR` | 2 new + NP 447-461 | 2.9e-03 | 3 | 3 | RSIR, RFMR, A149V | |
| `MRDNQQRASAGQISIQPTFSVQR` | 3 new + NP 397-416 | 4.4e-03 | 3 | 3 | RSIG, AFMR | |

#### Substitution (9)

| Peptide | Viral origin | Best PEP | Spectra | Runs | Conditions | Note |
|---|---|---|---|---|---|---|
| `VADQELGDAPFLDRR` | NS1 L36R | 9.2e-88 | 107 | 43 | RSIG, RSIR, RFMG, RFMR, AFMG, AFMR, A149V | see 6.3 |
| `REDGQQTPKLEHR` | PB1-F2 Q22E | 4.7e-03 | 8 | 8 | RSIG, AFMG, AFMR | ordinary PB1-F2 of the strain, see 6.3 |
| `TLDFHDSNVKK` | HA N460K | 6.8e-03 | 10 | 8 | RSIR, RFMG, RFMR, A149V | |
| `NAVGVLIGGLER` | NS1 W187R | 4.8e-05 | 8 | 7 | RSIG | |
| `SQLKNNAKEIGK` | HA N478K | 8.8e-03 | 6 | 6 | RFMG, RFMR, A149V | |
| `VSTRGVQIASNENMETMESSTLELR` | NP L358V | 2.6e-03 | 5 | 5 | RSIG, RSIR, RFMR | |
| `MLPFDRTTIMAAFNGNTEGR` | NP N417M | 9.1e-03 | 6 | 5 | RSIG, AFMG, AFMR, A149V | |
| `MSRNPGNAEFEDLTFLAR` | NP E244M | 2.3e-03 | 4 | 4 | RFMG, AFMR | |
| `GVGTMVMELVK` | NP R195K | 2.2e-03 | 3 | 3 | RFMR, AFMG, A149V | |

#### No viral anchor (4)

| Peptide | Viral origin | Best PEP | Spectra | Runs | Conditions | Note |
|---|---|---|---|---|---|---|
| `FHWEGLQDFISK` | longest viral match 3 of 12 residues | 1.5e-06 | 20 | 9 | RSIG, RSIR | same DPR as `NMAHWRVPQR` |
| `NMAHWRVPQR` | longest viral match 3 of 10 residues | 4.5e-04 | 14 | 7 | RSIG | same DPR as `FHWEGLQDFISK` |
| `NPALHSGQR` | longest viral match 4 of 9 residues | 4.7e-03 | 7 | 4 | RFMG | same DPR as `RQEIWASTK` |
| `RQEIWASTK` | longest viral match 3 of 9 residues | 5.6e-04 | 3 | 3 | RFMG | same DPR as `NPALHSGQR` |

![Figure 2: runs with an identification, per peptide and condition](figures/dpr_peptides_by_condition.png)

*Figure 2. Number of runs per condition in which each peptide was identified. The top block shows the 37 peptides that pass all filters. The lower blocks show peptides that fail exactly one filter, grouped by that filter.*

### 6.3 A closer look at the substitution peptides

This check goes beyond the automated analysis. For each substitution peptide, we counted how many DPR sequences carry the variant residue and how many carry the residue found in the reference PR8 protein (UniProt).

**`REDGQQTPKLEHR` (PB1-F2, E instead of Q22) is not DPR evidence.** 369 DPR sequences contain this E version of the peptide, and none contains the reference Q version. The PR8 genome from which the DPR sequences were translated thus encodes E at this position, unlike the PR8 reference protein in our database. The reference version of the peptide (`RQDGQQTPKLEHR`) is not detected at all in the data. The peptide is therefore most likely the ordinary PB1-F2 of the PR8 strain used, and it looks DPR-specific only because the reference database carries Q here. In addition, an encoded E cannot be told apart from a deamidated Q (Section 2.3). Filter 7 checks only the opposite direction, so this case slipped through. **Excluding it leaves 36 peptides.**

**For the other eight, the changed residue marks a junction.** The reference residue occurs in thousands of DPR sequences; NP `GVGTMVMELVR`, for instance, is found in 6,321. The variant occurs in only 1 to 23. The variant residue is therefore not a feature of the underlying genome. It appears in those DPR sequences where a junction falls exactly at this position, so that the first residue after the junction differs from the reference. Consistent with this, the changed residue is always the first or the last residue of the peptide.

These peptides are effectively junction peptides with only a single residue beyond the junction, which makes them weak evidence. A point mutation in a virus stock would produce exactly the same peptide, and each PR8 condition here is a separately generated mutant virus. Two pieces of information would settle the question: the DelVG junction coordinates, and the genome sequences of the virus stocks.

**`VADQELGDAPFLDRR` (NS1, R instead of L36) needs separate scrutiny.** It is by far the strongest identification (PEP ~10⁻⁸⁸, 107 spectra), found in every one of the 43 PR8 runs. Yet it occurs in only one DPR sequence. A DPR from a single DelVG would not be expected to appear so uniformly across seven different virus stocks. A sequence feature shared by all PR8 stocks, or a sample-processing artefact, seems at least as likely. It should not serve as the lead example without further checks.

### 6.4 Patterns across conditions

- **Most core peptides occur in several PR8 conditions.** For example, `MQLCVYGPAVASGYDFER` appears in 33 runs across all seven PR8 conditions, and `RMENFLDIWTYNAELLVLLENER` in 20 runs. By construction (filter 4), none occurs in mock.
- **Some peptides are restricted to one condition.** `GRWTTNTETEK` occurs only in A149V and `NAVGVLIGGLER` only in RSIG. The four peptides without a viral anchor form two pairs. Each pair comes from the same DPR sequence and was found exclusively or mostly in a single condition: `FHWEGLQDFISK` + `NMAHWRVPQR` in RSIG, and `NPALHSGQR` + `RQEIWASTK` in RFMG. Two independent peptides from one DPR in the same condition is a consistent pattern. DelVG populations differ between virus stocks, so condition-specific DPRs are plausible. A human protein missing from the database would instead be expected in mock and in the other conditions as well. These are the only two DPR sequences supported by more than one peptide.
- **Three peptides were also sequenced in H7N7 SC35M-infected cells:** `EEIRRENMETMESSTLELR`, `QTYDWTLNDNMTK` and `FYIQMCLRSR`. This dataset shows low-level carryover between runs (traces of PR8 peptides appear in mock samples), which could explain occasional hits. But the second half of `EEIRRENMETMESSTLELR` does not exist in SC35M NP, so this peptide cannot come from an SC35M DelVG. `QTYDWTLNDNMTK` was found in more SC35M runs (3) than in any PR8 condition. We regard these three as doubtful.

### 6.5 One peptide per DPR

A peptide is usually contained in many overlapping DPR sequences: in total, 1,128 DPR sequences contain at least one of the 37 peptides. The smallest set of DPR sequences that explains all 37 peptides has 35 members, and only the two named above contain more than one peptide. In proteomics, a protein identified by a single peptide is regarded as tentative, because the entire identification rests on one sequence.

## 7. How confident can we be?

- **Random matches.** One decoy remains among 37 peptides, so about one false hit is expected. Because the count is so small, the estimate is imprecise. Statistically, one observed decoy is compatible with up to about five false hits (95 % upper bound), that is, up to roughly 13 % of the list.
- **Systematic errors are not covered by the decoys.** Decoys mimic random matches only. They do not capture spectra that come from human peptides missing from the database, from chemically modified viral or human peptides, or from mixtures of two peptides, when such a spectrum happens to fit a DPR sequence. Filters 4, 6, 7 and 9 target these cases, but they are not exhaustive: Section 6.3 describes one case that slipped through.
- **Partial matches.** Many core peptides consist mostly of ordinary viral sequence and differ only in a short junction or new part. The ordinary part can dominate the spectrum. What matters is whether fragments covering the junction or the new residues are present. This has not been checked yet.
- **Not yet done:**
  - inspecting the spectra
  - comparing peptide junctions with the DelVG junctions observed by sequencing (the DPR database does not record the junction position of each sequence)
  - checking the four peptides without a viral anchor against newer human sequence databases and known variants.

## 8. Side effect: the ordinary viral proteins

MaxQuant does not report peptides directly. It assembles them into **protein groups**. A peptide shared by several proteins is assigned to the protein that explains the most peptides, and proteins without evidence of their own are merged into a group or dropped. Many DPR sequences contain almost a whole viral protein plus a few extra residues. As a result, they absorb the peptides of the ordinary viral proteins:

| Viral protein | This search (with DPRs, 67 runs) | Standard search (42 runs) |
|---|---|---|
| PB2 | group of 5,017 sequences, 1,260 MS/MS | group of 2, 3,179 MS/MS |
| PB1 | group of 1,345 sequences, 1,561 MS/MS | group of 1, 2,771 MS/MS |
| PB1-F2 | group of 591 sequences, 23 MS/MS | group of 1, 35 MS/MS |
| PA, PA-X | group of 3,659 sequences, 804 MS/MS | group of 2, 2,405 MS/MS |
| HA | group of 524 sequences, 2,264 MS/MS | group of 1, 1,957 MS/MS |
| NP | group of 338 sequences, 16,702 MS/MS | group of 1, 11,372 MS/MS |
| NA | group of 119 sequences, 45 MS/MS | group of 1, 39 MS/MS |
| M1 | **not reported** | group of 1, 3,752 MS/MS |
| M2 | group of 3 sequences, 7 MS/MS | group of 2, 951 MS/MS |
| NS1 | group of 325 sequences, 87 MS/MS | group of 1, 5,674 MS/MS |
| NEP | group of 3 sequences, 63 MS/MS | group of 2, 763 MS/MS |

The standard search has fewer runs, so its counts should be *lower*, not higher. Viral protein levels from this search are therefore not meaningful. **This search should be used for the DPR question only; viral protein quantification should come from the standard search.**

## 9. Proposed next steps

1. **Inspect the spectra** of the 24 core peptides, starting with those found in many runs. Check specifically for fragments that cover the junction or the new residues.
2. **Compare with DelVG sequencing.** Check whether each peptide's switch point coincides with a DelVG junction actually observed by sequencing. A peptide whose junction matches an observed DelVG is far more convincing. This requires the junction coordinates for each DPR sequence.
3. **Obtain the genome sequences of the virus stocks**, to distinguish point variants from junctions (Section 6.3).
4. **Validate the best candidates with synthetic peptides**, ideally isotope-labelled. Identical fragmentation and retention time is the accepted standard of proof. A targeted measurement (PRM) of these peptides would also be more sensitive than the present data.
5. **Refine the analysis script** so that it automatically flags Q→E and N→D substitutions (deamidated known peptides) and differences between the genome behind the DPR sequences and the PR8 reference proteins.

---

## Glossary

| Term | Meaning |
|---|---|
| Carryover | traces of one sample appearing in the next measurement on the same instrument |
| Chimeric spectrum | an MS/MS spectrum containing fragments of two peptides that were isolated together |
| Contaminants | proteins introduced during sample handling (keratins, trypsin, serum proteins); always included in the search so that they are not misassigned |
| Database | the list of protein sequences MaxQuant searches; only these can be identified |
| Deamidation | chemical change of N to D or Q to E during sample preparation; the product is identical to an encoded D or E |
| Decoy | a reversed protein sequence added to the search as a negative control; every decoy match is false |
| FDR (false discovery rate) | the expected share of false identifications in an accepted list; the standard is 1 % |
| Fragment / MS/MS spectrum | masses of the pieces of one peptide after it is broken apart; the "fingerprint" of its sequence |
| Isobaric | having the same mass (for example, G + V and R) |
| m/z | mass divided by charge; what a mass spectrometer measures |
| Main search | MaxQuant's regular search, as opposed to the second-peptide search |
| Match between runs | transfer of an identification from one run to another based on mass and retention time only, without a fragment spectrum |
| Missed cleavage | a K or R site that trypsin did not cut |
| PEP (posterior error probability) | estimated probability that one specific identification is wrong |
| ppm (parts per million) | unit of mass accuracy; 4.5 ppm of 1,500 Da is 0.007 Da |
| Protein group | proteins that MaxQuant cannot distinguish with the peptides found and therefore reports together |
| Retention time | the time at which a peptide leaves the chromatography column; reproducible for a given peptide |
| Run | one measurement of one sample |
| Second-peptide search | MaxQuant's search for a second peptide in a spectrum already explained by another |
| Target | a real protein sequence in the database, as opposed to a decoy |
| Tryptic peptide | a peptide produced by trypsin, ending in K or R |

## Files

- `answers.md` contains the complete, automatically generated results. In addition to the tables above, it lists the 60 peptides that fail exactly one filter and the peptides lost at each step.
- `figures/` and `tables/` contain the figures and per-peptide tables.
- `scripts/dpr_evidence.py` is the analysis script. Its header documents all definitions.
