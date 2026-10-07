# Pre-registration — clean rerun of the headline T-RIP RNA model, and DCP1's score — v01

Written 2026-10-07 18:05 UTC. No labels or expression values had been computed in this run when this was written. It follows HEADLINEAUDIT_v01 §4.

## 1. Data (fixed)
- **Labels.** These follow the DECAPLOCK_v01 §1 frozen rule.
  - S-BSST1096 IP libraries (DCP1/GFP × NS/HS × 2), complete files.
  - kallisto 0.51.1 against Ensembl Plants 63 TAIR10 cDNA + ncRNA; transcript est_counts are summed per gene.
  - CPM on n_processed. log2FC = log2((mean DCP1 CPM + 1)/(mean GFP CPM + 1)).
  - Universe: mean CPM ≥1 in both arms, minus ATCG/ATMG and the AT2G07590–07830 block.
  - **Gate G1:** n 5,900 / 4,824 (±1 %) and median +0.409 / +0.573 (±0.01). If G1 fails, the run stops and is reported.
- **Expression.** S-BSST1096 total RNA, lines 489 and 490, first 10 M read pairs, same index. E = log2(TPM + 1).
- **Features.** The 12 sequence features of TOTALLOCK (log10 lengths, exon count, GC, GC3, proximal-600-nt U/A/Urun/U5/CA), taken from `recovered/2026_ERCO2COREGULONS_DCP1FEATURES_REGEN_v01.tsv` (TAIR10 iGenomes, `.1` models). **Deviation:** the original used Ensembl 63 features.
- **Groups.** cd-hit-est −c 0.8 −n 5 on the longest Ensembl 63 cDNA per gene.

## 2. Changes from the headline (HEADLINEAUDIT_v01)
1. Expression from the **GFP line (489) only**: E489_cond and dE489 = E489_HS − E489_NS.
2. **No line-difference term.**
3. **DCP1 (AT1G08370) removed from the universe.** It is scored only by the fold models.

## 3. Model and evaluation (as TOTALLOCK §4)
- HistGradientBoostingRegressor (max_iter 300, lr 0.05, 15 leaves, l2 1.0) on the continuous log2FC.
- GroupKFold 5 × 5 repeats over clusters; out-of-fold predictions averaged.
- Primary metric: Spearman. Secondary: decile AUC (top vs bottom 10 % of y).
- Feature sets:
  - `seq`;
  - `clean` = seq + E489 + dE489;
  - `reproduce` = the original seq + Eboth + dEboth + L. This is a reproduction check only.
- **H-C1:** ΔSpearman(clean − seq) > 0, one-sided p < 0.05 by 2,000-sample cluster bootstrap, NS and HS separately (BH across the two).
- **Null:** 100 label permutations, clean feature set, decile AUC. The clean AUC must exceed the null 95th percentile.
- **Reproduction:** `reproduce` should fall within ±0.03 of 0.804 / 0.851. If it does not, that is reported as a deviation (the features differ, see §1).

## 4. DCP1 readout (descriptive)
For the clean model, report:
- predicted log2FC (mean of 25 fold models, with range);
- percentile among out-of-fold predictions;
- fraction of bottom-decile genes DCP1 outranks (the gene-level AUC equivalent) and fraction of top-decile genes it outranks;
- DCP1's observed log2FC and rank, and the G4c transgene excess.

Code: `2026_ERCO2COREGULONS_CLEANPIPE_v01.sh` (sha256 959260e4…25f8) and `2026_ERCO2COREGULONS_CLEANRUN_v01.py` (sha256 9e5e0c80…bd39). Not edited after the first label is computed; any bug fix is logged as a deviation.
