"""CLEANLOCK_v01 single run. Inputs (in ./work): ip_quant/<S>/abundance.tsv (8 IP libs, full),
tquant/<S>/abundance.tsv (11 total libs, first 10 M pairs), tx2gene.tsv, rep_c80.clstr,
features = 2026_ERCO2COREGULONS_DCP1FEATURES_REGEN_v01.tsv. Writes CLEANRUN_results.json + tables."""
import json, re, sys, numpy as np, pandas as pd, hashlib
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import roc_auc_score
W = 'work/'; BAIT = 'AT1G08370'; RNG = np.random.default_rng(20261007)
t2g = pd.read_csv(W + 'tx2gene.tsv', sep='\t', header=None, index_col=0)[1].to_dict()
def gene_table(d, col):
    out = {}
    for s in sorted(p for p in __import__('os').listdir(W + d) if not p.endswith('.log')):
        a = pd.read_csv(f'{W}{d}/{s}/abundance.tsv', sep='\t', index_col=0)
        out[s] = a[col].groupby(a.index.map(lambda t: t2g.get(t, t.split('.')[0]))).sum()
    T = pd.DataFrame(out).fillna(0); T.index = T.index.astype(str).str.upper(); return T
OUT = {}
# ---------- labels (DECAPLOCK §1 frozen rule) ----------
C = gene_table('ip_quant', 'est_counts')        # columns DCP1_NS_1 ... GFP_HS_2
proc = {}
for s in C.columns:
    js = json.load(open(f'{W}ip_quant/{s}/run_info.json')); proc[s] = js['n_processed']
CPM = C.div(pd.Series(proc), axis=1) * 1e6
LAB = {}
for cond in ('NS', 'HS'):
    d = [c for c in CPM if c.startswith('DCP1_' + cond)]; g = [c for c in CPM if c.startswith('GFP_' + cond)]
    D, G = CPM[d].mean(1), CPM[g].mean(1); lfc = np.log2((D + 1) / (G + 1))
    num = lambda x: int(x[2:7]) if x.startswith('AT2G') else -1
    ok = (D >= 1) & (G >= 1) & ~lfc.index.str.startswith(('ATCG', 'ATMG')) & ~lfc.index.map(lambda x: 7590 <= num(x) <= 7830)
    u = lfc[ok]; LAB[cond] = u
    OUT[f'G1_{cond}'] = dict(n=int(len(u)), median=round(float(u.median()), 4),
                            bait_lfc=round(float(lfc.get(BAIT, np.nan)), 3), bait_rank=int((u > u.get(BAIT, np.inf)).sum()) + 1)
anch = {'NS': (5900, 0.409), 'HS': (4824, 0.573)}
G1 = all(abs(OUT[f'G1_{c}']['n'] - anch[c][0]) <= 0.01 * anch[c][0] and abs(OUT[f'G1_{c}']['median'] - anch[c][1]) <= 0.01 for c in anch)
OUT['G1_pass'] = bool(G1); print('G1', OUT['G1_NS'], OUT['G1_HS'], 'PASS' if G1 else 'FAIL', flush=True)
pd.DataFrame(LAB).to_csv('2026_ERCO2COREGULONS_CLEANLABELS_v01.tsv', sep='\t')
# ---------- expression (TOTALLOCK §1) ----------
T = gene_table('tquant', 'tpm'); E = np.log2(T + 1); E.to_csv('2026_ERCO2COREGULONS_CLEANTOTALLOGTPM_v01.tsv', sep='\t')
cols = lambda line, cond: [c for c in E if c.startswith(f'{line}_{cond}')]
OUT['G4c_DCP1_490minus489'] = {c: round(float(E.loc[BAIT, cols('490', c)].mean() - E.loc[BAIT, cols('489', c)].mean()), 2) for c in ('NS', 'HS')}
E489 = {c: E[cols('489', c)].mean(1) for c in ('NS', 'HS')}
Eboth = {c: E[cols('489', c) + cols('490', c)].mean(1) for c in ('NS', 'HS')}
Lline = {c: E[cols('490', c)].mean(1) - E[cols('489', c)].mean(1) for c in ('NS', 'HS')}
# ---------- features + groups ----------
F = pd.read_csv('2026_ERCO2COREGULONS_DCP1FEATURES_REGEN_v01.tsv', sep='\t', index_col=0)
for c in ['L_tot', 'L_cds', 'L_5', 'L_3']: F[c] = np.log10(F[c] + 1)
SEQ = ['L_tot', 'L_cds', 'L_5', 'L_3', 'n_exon', 'GC', 'GC3', 'U3', 'A3', 'Urun', 'U5', 'CA']
clu = {}; cid = -1
for ln in open(W + 'rep_c80.clstr'):
    if ln.startswith('>Cluster'): cid += 1
    else: clu[ln.split('>')[1].split('...')[0].upper()] = cid
def fit_oof(X, y, groups, xb, reps=5, seed0=0):
    P = np.zeros(len(y)); B = []
    for r in range(reps):
        ug = np.unique(groups); perm = dict(zip(ug, np.random.default_rng(seed0 + r).permutation(len(ug))))
        g2 = np.array([perm[g] for g in groups])
        for tr, te in GroupKFold(5).split(X, y, g2):
            m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, l2_regularization=1.0, random_state=r)
            m.fit(X[tr], y[tr]); P[te] += m.predict(X[te])
            if xb is not None: B.append(float(m.predict(xb)[0]))
    return P / reps, B
def dauc(s, y, q=0.1):
    lo, hi = np.quantile(y, [q, 1 - q]); m = (y <= lo) | (y >= hi)
    return roc_auc_score((y[m] >= hi).astype(int), s[m])
for cond in ('NS', 'HS'):
    u = LAB[cond]; u = u[u.index.isin(F.index) & u.index.isin(list(clu))]
    bait_in = BAIT in u.index
    u = u.drop(BAIT, errors='ignore')                    # bait excluded from universe (CLEANLOCK §2)
    X = F.loc[u.index, SEQ].copy()
    X['E489'] = E489[cond].reindex(u.index); X['dE489'] = (E489['HS'] - E489['NS']).reindex(u.index)
    X['Eboth'] = Eboth[cond].reindex(u.index); X['dEboth'] = (Eboth['HS'] - Eboth['NS']).reindex(u.index); X['Lline'] = Lline[cond].reindex(u.index)
    xb_all = F.loc[[BAIT], SEQ].copy()
    xb_all['E489'] = E489[cond][BAIT]; xb_all['dE489'] = E489['HS'][BAIT] - E489['NS'][BAIT]
    xb_all['Eboth'] = Eboth[cond][BAIT]; xb_all['dEboth'] = Eboth['HS'][BAIT] - Eboth['NS'][BAIT]; xb_all['Lline'] = Lline[cond][BAIT]
    y = u.values; grp = np.array([clu[g] for g in u.index]); R = dict(n=int(len(y)), bait_in_universe=bait_in)
    SETS = {'seq': SEQ, 'clean_seq+E489': SEQ + ['E489', 'dE489'], 'reproduce_seq+Eboth+L': SEQ + ['Eboth', 'dEboth', 'Lline']}
    preds = {}
    for name, fs in SETS.items():
        P, B = fit_oof(X[fs].values, y, grp, xb_all[fs].values)
        preds[name] = P
        R[name] = dict(spearman=round(spearmanr(P, y).correlation, 4), decile_auc=round(dauc(P, y), 4))
        if name == 'clean_seq+E489':
            b = float(np.mean(B)); lo, hi = np.quantile(y, [.1, .9])
            neg = P[y <= lo]; pos = P[y >= hi]
            R['DCP1'] = dict(pred_lfc=round(b, 3), pred_range=[round(min(B), 3), round(max(B), 3)],
                pct_among_oof=round(float((P < b).mean() * 100), 1), outranks_depleted=round(float((neg < b).mean()), 3),
                outranks_enriched=round(float((pos < b).mean()), 3), features={k: round(float(v), 4) for k, v in xb_all[fs].iloc[0].items()})
            np.save(f'oof_clean_{cond}.npy', P)
    # cluster bootstrap CIs for clean model and delta vs seq
    ug = np.unique(grp); gi = {q: np.where(grp == q)[0] for q in ug}; bs = []
    for _ in range(2000):
        idx = np.concatenate([gi[q] for q in RNG.choice(ug, len(ug))])
        bs.append((spearmanr(preds['clean_seq+E489'][idx], y[idx]).correlation - spearmanr(preds['seq'][idx], y[idx]).correlation,
                   dauc(preds['clean_seq+E489'][idx], y[idx])))
    bs = np.array(bs); R['clean_delta_spearman_ci'] = [round(float(np.percentile(bs[:, 0], q)), 4) for q in (2.5, 97.5)]
    R['clean_delta_p_one_sided'] = round(float((bs[:, 0] <= 0).mean()), 4); R['clean_decile_auc_ci'] = [round(float(np.percentile(bs[:, 1], q)), 4) for q in (2.5, 97.5)]
    # label-shuffle null (100 permutations, 1 repeat each, clean feature set)
    nul = []
    for p in range(100):
        ys = np.random.default_rng(1000 + p).permutation(y)
        Pn, _ = fit_oof(X[SETS['clean_seq+E489']].values, ys, grp, None, reps=1, seed0=1000 + p)
        nul.append(dauc(Pn, ys))
    R['null_decile_auc'] = dict(mean=round(float(np.mean(nul)), 4), q95=round(float(np.quantile(nul, .95)), 4), max=round(float(np.max(nul)), 4))
    OUT[cond] = R; print(cond, json.dumps(R), flush=True)
json.dump(OUT, open('2026_ERCO2COREGULONS_CLEANRUN_results_v01.json', 'w'), indent=1)
print('CLEANRUN DONE', flush=True)
