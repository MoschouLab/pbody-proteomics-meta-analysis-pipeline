# Clean rerun kit (CLEANLOCK_v01)

Pre-registered rerun of the headline T-RIP RNA model with GFP-line (489) expression only, no line-difference term,
DCP1 (AT1G08370) out of the universe. The lock is `2026_ERCO2COREGULONS_CLEANLOCK_v01.md`; do not edit the two locked
scripts (sha256 959260e4…25f8 and 9e5e0c80…bd39). Any bug fix is logged as a deviation.

Needs EBI access (ftp.ebi.ac.uk, www.ebi.ac.uk, ftp.ensemblgenomes.ebi.ac.uk), kallisto 0.51.1, cd-hit-est,
Python 3 with numpy, pandas, scipy, scikit-learn.

1. `bash 2026_ERCO2COREGULONS_CLEANPIPE_v01.sh` builds work/ref, tx2gene.tsv, rep_c80.clstr, study.json.
2. IP libraries (8, complete files) into `work/ip_quant/{DCP1,GFP}_{NS,HS}_{1,2}` — recipe in
   `2026_ERCO2COREGULONS_COUNTS_PIPE_v01.sh` (paths adapted to work/, index work/ref/txome.idx).
3. Total RNA (11 libraries of lines 489/490) into `work/tquant/<line>_<cond>_<n>`, first 10 M read pairs
   (40,000,000 FASTQ lines per mate) — recipe in `run_total.sh`, library list from `mk_total_libs.py`.
4. `python3 2026_ERCO2COREGULONS_CLEANRUN_v01.py` with `2026_ERCO2COREGULONS_DCP1FEATURES_REGEN_v01.tsv` in the
   same folder. Stop if gate G1 fails (n 5,900 / 4,824 ±1 %, median log2FC +0.409 / +0.573 ±0.01).
