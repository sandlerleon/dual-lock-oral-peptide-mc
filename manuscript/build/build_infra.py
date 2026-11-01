# -*- coding: utf-8 -*-
"""Infrastructure for build_manuscript.py: author-year citations (Springer basic style), docx building with the shared engine, number formatting."""
import json
import os
import re

import docx_engine as E


class Harvard:
    """Elsevier author-year style. [@a;@b] -> (Author et al., 2016; Author, 2018);  [!@a] -> Author et al. (2016)."""
    PAT = re.compile(r"\[(!?)@([^\]]+)\]")

    def __init__(self, lib):
        self.lib, self.used = lib, set()

    def _keys(self, m):
        return [x.strip().lstrip("@") for x in m.group(2).split(";")]

    def scan(self, text):
        for m in self.PAT.finditer(text):
            for k in self._keys(m):
                if k not in self.lib:
                    raise KeyError(k)
                self.used.add(k)

    def resolve(self, text):
        def rep(m):
            keys = self._keys(m)
            if m.group(1) == "!":
                k = keys[0]; return "%s (%d)" % (self.lib[k][0], self.lib[k][1])
            return "(" + "; ".join("%s, %d" % (self.lib[k][0], self.lib[k][1]) for k in keys) + ")"
        return self.PAT.sub(rep, text)

    def entries(self):
        unused = [k for k in self.lib if k not in self.used]
        if unused:
            raise ValueError("never cited: " + ", ".join(unused))
        return sorted((self.lib[k][2], self.lib[k][3]) for k in self.used)


def n(x, d=2):
    return ("{:.%df}" % d).format(x).replace("-", "−")


def sgn(x, d=2):
    return ("{:+.%df}" % d).format(x).replace("-", "−")


def ci(t, d=2):
    """[value, lo, hi] -> 'v (lo to hi)' """
    return "%s (%s to %s)" % (n(t[0], d), n(t[1], d), n(t[2], d))


def md_inline(text):
    """**bold** and *italic* markers inside block text -> engine runs. The engine understands $math$ and [@cite] only, so bold/italic are
    handled by build_manuscript via split()."""
    return text


def build(path, blocks, cit, refs=True):
    doc = E.Doc(cit)
    for b in blocks:
        k = b[0]
        if k == "title":
            doc.title(b[1])
        elif k == "centered":
            doc.centered(b[1], b[2] if len(b) > 2 else 22)
        elif k == "blank":
            doc.blank()
        elif k == "h1":
            doc.h1(b[1])
        elif k == "h2":
            doc.h2(b[1])
        elif k == "p":
            doc.p(b[1])
        elif k == "p0":
            doc.p(b[1], indent=False)
        elif k == "bullet":
            doc.p("• " + b[1], indent=False)
        elif k == "eq":
            doc.eq(b[1], b[2])
        elif k == "fig":
            doc.fig(b[1], b[2], b[3])
        elif k == "table":
            doc.table(b[1], b[2], b[3], b[4] if len(b) > 4 else None, size=b[5] if len(b) > 5 else 16)
        elif k == "refs":
            doc.h1("References")
            for _, text in cit.entries():
                par = doc.d.add_paragraph()
                E._spacing(par, 0, 60)
                par.paragraph_format.left_indent = E.Emu(317 * 635)
                par.paragraph_format.first_line_indent = E.Emu(-317 * 635)
                E.add_runs(par, text, 18)
    doc.save(path)


def scan_blocks(cit, blocks):
    for b in blocks:
        for x in b[1:]:
            if isinstance(x, str):
                cit.scan(x)
            elif isinstance(x, list):
                for row in x:
                    for c in (row if isinstance(row, list) else [row]):
                        if isinstance(c, str):
                            cit.scan(c)
