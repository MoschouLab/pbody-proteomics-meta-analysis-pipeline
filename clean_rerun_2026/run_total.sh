#!/bin/bash
B=https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/096/S-BSST1096/Files
cd /home/claude/rebuild; mkdir -p tquant
while IFS=$'\t' read s f1 f2; do
  [ -f tquant/$s/abundance.tsv ] && { echo "$s done"; continue; }
  rm -f /tmp/r1_$s /tmp/r2_$s; mkfifo /tmp/r1_$s /tmp/r2_$s
  (wget -q -O - "$B/$f1" | funzip 2>/dev/null | head -n 40000000 > /tmp/r1_$s) &
  (wget -q -O - "$B/$f2" | funzip 2>/dev/null | head -n 40000000 > /tmp/r2_$s) &
  ./kallisto/kallisto quant -i tair10_e63.idx -o tquant/$s -t 2 /tmp/r1_$s /tmp/r2_$s > tquant/$s.log 2>&1
  echo "$s exit $? $(date +%T) $(grep -E 'processed' tquant/$s.log | tail -1)"
  wait; rm -f /tmp/r1_$s /tmp/r2_$s
done < total_libs.tsv
echo ALLDONE
