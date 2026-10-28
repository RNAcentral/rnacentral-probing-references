#!/usr/bin/env python3
# Small-RNA test transcriptome (tRNAs + rDNA unit, 5S, mt 16S); from the repo root: build_trna_rrna_test.py 28
import gzip
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

# the two rDNA units the rf-eval windows sit on, 5S, mt 16S and mt 12S; standalone 18S/5.8S/28S would duplicate them
RRNA = ['URS0000ABD879_9606', 'URS0000ABD87F_9606', 'URS00000F9D45_9606', 'URS000080E357_9606', 'URS000044DFF6_9606']
STRUCTURES = ['rfeval/rrna_common.db', '../../scripts/homo_sapiens/rnastruct00040_trna.db']

release = sys.argv[1]
org = Path(f'release{release}/homo_sapiens')
work = org / 'work'
work.mkdir(exist_ok=True)

trna = set()
for line in gzip.open(org / 'reference/homo_sapiens.GRCh38.filtered.gtf.gz', 'rt'):
    f = line.split('\t')
    if f[2] == 'transcript' and 'transcript_biotype "tRNA"' in f[8]:
        trna.add(re.search(r'transcript_id "([^"]+)"', f[8]).group(1).rsplit('.', 1)[0])

cache_path = work / 'trna_rrna_seqs.json'
cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
# GtRNAdb's export covers most tRNAs in one download; the API fills in the rest (mt tRNAs, rRNA)
gtrnadb = work / 'gtrnadb.fasta.gz'
if not gtrnadb.exists():
    urllib.request.urlretrieve(f'https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/{release}.0/sequences/by-database/gtrnadb.fasta.gz', gtrnadb)
for rec in gzip.open(gtrnadb, 'rt').read().split('>')[1:]:
    head, *body = rec.split('\n')
    if head.split()[0] in trna:
        cache[head.split()[0]] = ''.join(body).upper().replace('U', 'T')
for urs in sorted(trna) + RRNA:
    if urs in cache:
        continue
    for wait in (1, 5, 20, 60):
        time.sleep(wait)
        req = urllib.request.Request(f'https://rnacentral.org/api/v1/rna/{urs.split("_")[0]}?format=json')
        try:
            body = urllib.request.urlopen(req, timeout=60).read().decode()
            cache[urs] = json.loads(body)['sequence'].upper().replace('U', 'T')
            break
        except (ValueError, OSError):
            continue  # throttled responses come back as the HTML home page
    else:
        sys.exit(f'{urs}: no sequence from rnacentral.org')
    cache_path.write_text(json.dumps(cache, indent=0, sort_keys=True))

seqs = {u: cache[u] for u in sorted(trna)}
# a record holding another whole (e.g. a pre-tRNA with leader) would split that tRNA's reads
contained = {a for a in seqs for b in seqs if a != b and seqs[b] in seqs[a]}
records = [(u, s) for u, s in seqs.items() if u not in contained] + [(u, cache[u]) for u in RRNA]
assert all(not re.search('[^ACGTN]', s) for _, s in records)

fa = org / 'reference/homo_sapiens.trna_rrna.transcriptome.fa'
with open(fa, 'w') as fh:
    for rid, s in records:
        fh.write(f'>{rid}\n' + '\n'.join(s[i:i + 60] for i in range(0, len(s), 60)) + '\n')
print(f'{fa}: {len(records)} records ({len(trna)} tRNA URS, {len(contained)} dropped as containing another, {len(RRNA)} rRNA)')

# existing structures re-keyed onto the new records; Rfam structures are never rebuilt
db, win = [], ['# ref_seq_id\tstart\tend\tstructure_id']
for src in STRUCTURES:
    lines = (org / src).read_text().split('\n')
    for i in range(0, len(lines) - 2, 3):
        name, s, ss = lines[i][1:], lines[i + 1].replace('U', 'T'), lines[i + 2]
        hits = [(rid, r.find(s)) for rid, r in records if s in r]
        if not hits:
            print(f'  {name}: sequence not in any record, skipped')
            continue
        rid, pos = hits[0]
        sid = f'{rid}__{name.split("__", 1)[1]}'
        db.append(f'>{sid}\n{lines[i + 1]}\n{ss}')
        win.append(f'{rid}\t{pos + 1}\t{pos + len(s)}\t{sid}')
        print(f'  {name} -> {sid} {pos + 1}-{pos + len(s)}' + (f' ({len(hits)} records match)' if len(hits) > 1 else ''))
(org / 'rfeval/trna_rrna.db').write_text('\n'.join(db) + '\n')
(org / 'rfeval/trna_rrna.windows.tsv').write_text('\n'.join(win) + '\n')
