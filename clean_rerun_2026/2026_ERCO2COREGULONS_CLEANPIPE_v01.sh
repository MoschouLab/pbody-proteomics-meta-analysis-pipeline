#!/bin/bash
# CLEANLOCK_v01 data pipeline: Ensembl Plants 63 TAIR10 index, IP libs (full), total libs (first 10M pairs),
# rep cDNA clusters (cd-hit-est 0.8). Run from the CLEANRUN folder inside the conda env 'cleanrun'.
set -u; mkdir -p work/ip_quant work/tquant work/ref; cd work
E63=https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/arabidopsis_thaliana
BS=https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/096/S-BSST1096/Files
[ -s ref/cdna.fa.gz ] || curl -sfL -o ref/cdna.fa.gz $E63/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
[ -s ref/ncrna.fa.gz ] || curl -sfL -o ref/ncrna.fa.gz $E63/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz
[ -s ref/tx.fa.gz ] || cat ref/cdna.fa.gz ref/ncrna.fa.gz > ref/tx.fa.gz
[ -s ref/txome.idx ] || kallisto index -i ref/txome.idx ref/tx.fa.gz
python3 - <<'PY'
import gzip,re
from collections import defaultdict
L=defaultdict(lambda:("",0));seq={};k=None;b=[]
out=open("tx2gene.tsv","w")
with gzip.open("ref/tx.fa.gz","rt") as f:
    for ln in f:
        if ln[0]==">":
            if k: seq[k]="".join(b)
            k=ln[1:].split()[0];b=[];m=re.search(r"gene:(\S+)",ln);g=(m.group(1) if m else k.split(".")[0]).upper()
            out.write(f"{k}\t{g}\n");gene_of=g;seq.setdefault("__g__"+k,g)
        else: b.append(ln.strip())
    seq[k]="".join(b)
rep={}
for t,s in seq.items():
    if t.startswith("__g__"): continue
    g=seq["__g__"+t]
    if g.startswith("AT") and len(s)>len(rep.get(g,"")): rep[g]=s
with open("rep_cdna.fa","w") as f:
    for g,s in rep.items(): f.write(f">{g}\n{s}\n")
PY
[ -s rep_c80.clstr ] || { cd-hit-est -i rep_cdna.fa -o rep_c80 -c 0.8 -n 5 -M 8000 -T 4 -d 0 > cdhit.log 2>&1; }
curl -sfL -o study.json "https://www.ebi.ac.uk/biostudies/api/v1/studies/S-BSST1096/info" || true
curl -sfL -o filelist.json "https://www.ebi.ac.uk/biostudies/files/S-BSST1096/S-BSST1096.json" || true
echo "REF READY $(date +%T)"
