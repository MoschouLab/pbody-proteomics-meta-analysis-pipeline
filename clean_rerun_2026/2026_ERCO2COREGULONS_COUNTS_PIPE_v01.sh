#!/bin/bash
# Re-quantify the deposited T-RIP FASTQs (BioStudies S-BSST1096) into RAW COUNTS.
# One sample at a time: download -> kallisto quant -> delete FASTQs. Peak disk ~5 GB.
set -u
cd /home/claude
BASE="https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/096/S-BSST1096/Files"
IDX=ref/txome.idx
mkdir -p quant fq
python3 - <<'PY' > sample_map.tsv
import json
f=json.load(open("trip_files.json"))
by={}
for x in f:
    d=x["desc"]                       # e.g. "GFP NS 1 FORWARD"
    parts=d.split()
    key="_".join(parts[:3])           # GFP_NS_1
    by.setdefault(key,{})[parts[3]]=x["path"]
for k,v in sorted(by.items()):
    print(f"{k}\t{v['FORWARD']}\t{v['REVERSE']}")
PY
while IFS=$'\t' read -r S R1 R2; do
  if [ -f "quant/$S/abundance.tsv" ]; then echo "[skip] $S already quantified"; continue; fi
  echo "[$(date +%H:%M:%S)] $S  downloading"
  curl -sL -m 3600 -o "fq/${S}_1.fq.gz" "$BASE/$R1" || { echo "[fail dl] $S"; continue; }
  curl -sL -m 3600 -o "fq/${S}_2.fq.gz" "$BASE/$R2" || { echo "[fail dl] $S"; continue; }
  echo "[$(date +%H:%M:%S)] $S  downloaded $(du -sh fq | cut -f1); quantifying"
  kallisto quant -i "$IDX" -o "quant/$S" -t 2 "fq/${S}_1.fq.gz" "fq/${S}_2.fq.gz" > "quant/${S}.log" 2>&1
  if [ -f "quant/$S/abundance.tsv" ]; then
     echo "[$(date +%H:%M:%S)] $S  DONE  $(grep -E 'pseudoaligned' quant/${S}.log | head -1)"
  else
     echo "[$(date +%H:%M:%S)] $S  QUANT FAILED"; tail -3 "quant/${S}.log"
  fi
  rm -f "fq/${S}_1.fq.gz" "fq/${S}_2.fq.gz"
done < sample_map.tsv
echo "[$(date +%H:%M:%S)] pipeline finished"
