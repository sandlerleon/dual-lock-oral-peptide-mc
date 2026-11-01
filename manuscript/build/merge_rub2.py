# -*- coding: utf-8 -*-
"""Selective merge of the Rubriq (AJE) language edit for the Dual-Lock manuscript: keep the edited file (equations, tables, figures intact) and put back my paragraphs wherever the edit
changed the meaning or damaged the text (math symbols dropped, wrong subjects, garbled sentences).

    python merge_rub2.py <mine.docx> <rub.docx> <out.docx>
"""
import copy
import re
import sys

import docx
from docx.oxml.ns import qn

MINE, RUB, OUT = sys.argv[1:4]
KEYS = [
    "With the reference calibrated to 1.2%, Dual-Lock reached",
    "A comparison of SNAC with sodium caprate showed",
    "D-α-tocopheryl polyethylene glycol 1000 succinate (TPGS) inhibits",
    "With tclear",
    "Optimized design out of sample.",
    "With the reference arm calibrated to a median of 1.2%, the Dual-Lock arm P1",
    "Because all arms share the same draws",
    "Removing one component from P1",
    "Two results differ from the original design intuition",
    "With the assumed correlation set",
    "Fig. 5. Monte Carlo convergence",
    "Fig. 6. Benefit-risk analysis",
    "so that the peptide released into a small niche",
    "The qualitative ordering is the same",
    "Four comparisons with published data",
    "Limitations. All common priors",
    "Implications. The most useful outputs",
    "A mechanistic, evidence-graded Monte Carlo analysis of a co-localized",
    "All code, tests, numerical results (JSON and CSV)",
]
SPELL = [("co-localiz", "colocaliz"), ("Co-localiz", "Colocaliz"), ("modelled", "modeled"), ("artefact", "artifact"), ("wash-out", "washout"), ("damage-hours", "damage hours"), ("co-formulation", "coformulation")]


def nonempty(d):
    return [p for p in d.paragraphs if p.text.strip()]


def fix_spelling(par):
    for t in par._p.iter(qn("w:t")):
        if t.text:
            x = t.text
            for a, b in SPELL:
                x = x.replace(a, b)
            t.text = x


mine = docx.Document(MINE); rub = docx.Document(RUB)
pm, pr = nonempty(mine), nonempty(rub)
assert len(pm) == len(pr), (len(pm), len(pr))
done = []
for key in KEYS:
    hits = [i for i, p in enumerate(pm) if key in p.text]
    if not hits:
        print("NOT FOUND:", key); continue
    for i in hits:
        new = copy.deepcopy(pm[i]._p)
        pr[i]._p.addnext(new); pr[i]._p.getparent().remove(pr[i]._p)
        fix_spelling(docx.text.paragraph.Paragraph(new, pr[i]._parent))
        done.append(i)
print("restored paragraphs:", len(done))
rub.save(OUT)
