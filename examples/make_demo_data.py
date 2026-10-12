#!/usr/bin/env python3
"""
make_demo_data.py -- small SIMULATED dataset to try the plotting scripts (no real data).

Creates in examples/demo_data/:
  genome.fa (+.fai)       one 6 kb chromosome "chr1"
  genes.gtf               FAST (+ strand, high labeling) and STABLE (- strand, low labeling)
  pre.bam (+.bai)         pre-consensus reads with CB/UB tags, 2 cells
  consensus.bam (+.bai)   one consensus read per UMI (dynast rule, RN tag)
  mask_0based.csv         SNP mask (0-based), incl. a homozygous T>C SNP in FAST
"""
import os
import random
from collections import defaultdict

import pysam

random.seed(1)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "demo_data")
os.makedirs(out, exist_ok=True)

L = 6000
seq = list("".join(random.choice("ACGT") for _ in range(L)))
SNP = 2500                      # 0-based position of a genomic T>C SNP in FAST
seq[SNP] = "T"
seq = "".join(seq)
with open(f"{out}/genome.fa", "w") as f:
    f.write(">chr1\n" + "\n".join(seq[i:i + 60] for i in range(0, L, 60)) + "\n")
pysam.faidx(f"{out}/genome.fa")

with open(f"{out}/genes.gtf", "w") as f:
    for gid, name, s, e, st in (("G1", "FAST", 1001, 3000, "+"), ("G2", "STABLE", 3501, 5500, "-")):
        attr = f'gene_id "{gid}"; transcript_id "{gid}.1"; gene_name "{name}"; tag "Ensembl_canonical";'
        f.write(f"chr1\tdemo\tgene\t{s}\t{e}\t.\t{st}\t.\t{attr}\n")
        f.write(f"chr1\tdemo\ttranscript\t{s}\t{e}\t.\t{st}\t.\t{attr}\n")
        f.write(f"chr1\tdemo\texon\t{s}\t{e}\t.\t{st}\t.\t{attr}\n")

masked = {(SNP, "TC")}
for p in random.sample(range(1000, 5500), 30):
    if seq[p] in "TA":
        masked.add((p, "TC" if seq[p] == "T" else "AG"))
with open(f"{out}/mask_0based.csv", "w") as f:
    f.write("contig,genome_i,conversion\n")
    for p, c in sorted(masked):
        f.write(f"chr1,{p},{c}\n")

hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": L}]}
pre, umis = [], defaultdict(list)
k = 0
for cell, n_fast, n_stable in (("CELLAAAAAAAAAAAA", 40, 40), ("CELLCCCCCCCCCCCC", 12, 50)):
    for gstart, glen, rev, plab, n in ((2200, 500, False, 0.6, n_fast), (3600, 1800, True, 0.05, n_stable)):
        U, C = ("A", "G") if rev else ("T", "C")
        for u in range(n):
            umi = f"{'S' if rev else 'F'}{u:03d}{cell[4:8]}"
            st = gstart + random.randint(0, glen - 130)
            ts = [p for p in range(st, st + 130) if seq[p] == U]
            conv = set(random.sample(ts, min(len(ts), random.randint(1, 3)))) \
                if random.random() < plab else set()
            for r in range(random.choice([1, 1, 1, 2, 3, 5])):
                o = random.choice([0, 0, 10, 25])
                s0 = st + o
                s = list(seq[s0:s0 + 90])
                q = [37] * 90
                for i in range(90):
                    p = s0 + i
                    if p in conv or (p == SNP and not rev):
                        s[i] = C
                    elif seq[p] == U and random.random() < 0.0005:    # sequencing error
                        s[i] = C
                a = pysam.AlignedSegment()
                a.query_name = f"r{k}"; k += 1
                a.query_sequence = "".join(s)
                a.flag = 16 if rev else 0
                a.reference_id, a.reference_start, a.mapping_quality = 0, s0, 255
                a.cigar = [(0, 90)]
                a.query_qualities = pysam.qualitystring_to_array("".join(chr(x + 33) for x in q))
                a.set_tags([("CB", cell), ("UB", umi)])
                pre.append(a)
                umis[(cell, umi, rev)].append(a)


def write(path, reads):
    with pysam.AlignmentFile(path, "wb", header=hdr) as o:
        for a in sorted(reads, key=lambda r: r.reference_start):
            o.write(a)
    pysam.index(path)


write(f"{out}/pre.bam", pre)

cons = []                                   # dynast-style consensus: quality-weighted vote
for (cell, umi, rev), rs in umis.items():
    lo, hi = min(r.reference_start for r in rs), max(r.reference_end for r in rs)
    bs, qs = [], []
    for p in range(lo, hi):
        v = defaultdict(int)
        for r in rs:
            if r.reference_start <= p < r.reference_end:
                i = p - r.reference_start
                v[r.query_sequence[i]] += r.query_qualities[i]
        best = max(v.values()) if v else 0
        win = [b for b in v if v[b] == best]
        b = win[0] if v and best >= 27 and len(win) == 1 else seq[p]
        bs.append(b); qs.append(min(v.get(b, 0), 42))
    c = pysam.AlignedSegment()
    c.query_name = f"c_{umi}"; c.query_sequence = "".join(bs)
    c.flag = 16 if rev else 0
    c.reference_id, c.reference_start, c.mapping_quality = 0, lo, 255
    c.cigar = [(0, hi - lo)]
    c.query_qualities = pysam.qualitystring_to_array("".join(chr(x + 33) for x in qs))
    c.set_tags([("CB", cell), ("UB", umi), ("RN", len(rs))])
    cons.append(c)
write(f"{out}/consensus.bam", cons)
print(f"demo data written to {out}: {len(pre)} reads, {len(cons)} molecules, {len(masked)} masked sites")
