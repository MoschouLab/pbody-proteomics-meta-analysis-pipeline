import json, re
d = json.load(open("study.json")); files = {}
def walk(o):
    if isinstance(o, dict):
        p = o.get("path", "")
        if re.match(r"^(489|490)_", p): files[p] = o.get("size")
        for v in o.values(): walk(v)
    elif isinstance(o, list):
        for v in o: walk(v)
walk(d)
S = {}
for p in files:
    m = re.match(r"^(489|490)_(NS|HS)_(\d)_([12])", p); S.setdefault(m.group(1, 2, 3), {})[m.group(4)] = p
with open("total_libs.tsv", "w") as f:
    for k in sorted(S): f.write("\t".join(["_".join(k), S[k]["1"], S[k]["2"]]) + "\n")
print(open("total_libs.tsv").read())
