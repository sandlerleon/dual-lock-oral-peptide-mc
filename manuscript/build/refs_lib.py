# -*- coding: utf-8 -*-
"""Reference library (Elsevier author-year style). Metadata from Crossref (fetch_refs.py) checked on 1 Nov 2026; abstracts of the empirical sources were read in Europe PMC."""
import json, os
META = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "refs_meta.json"), encoding="utf-8"))
ABBR = {"The Journal of Pharmacology and Experimental Therapeutics": "J. Pharmacol. Exp. Ther.", "Science Translational Medicine": "Sci. Transl. Med.",
        "Journal of Controlled Release": "J. Control. Release", "Pharmaceutical Research": "Pharm. Res.", "International Journal of Pharmaceutics": "Int. J. Pharm.",
        "Advanced Drug Delivery Reviews": "Adv. Drug Deliv. Rev.", "Nanoscale Advances": "Nanoscale Adv.", "Pharmaceutics": "Pharmaceutics",
        "Clinical Pharmacokinetics": "Clin. Pharmacokinet.", "Scientific Reports": "Sci. Rep.", "Nature Reviews Drug Discovery": "Nat. Rev. Drug Discov.",
        "American Journal of Physiology-Gastrointestinal and Liver Physiology": "Am. J. Physiol. Gastrointest. Liver Physiol.", "Computer Physics Communications": "Comput. Phys. Commun.",
        "Computer Physics Communication": "Comput. Phys. Commun.", "Mathematics and Computers in Simulation": "Math. Comput. Simul.",
        "Communications in Statistics - Simulation and Computation": "Commun. Stat. Simul. Comput.", "Genetics": "Genetics", "Journal of Theoretical Biology": "J. Theor. Biol."}
YEAR = {"Granhall2019": 2019, "Drucker2020": 2020}


TITLE_OVERRIDE = {"Dintaman1999": "Inhibition of P-glycoprotein by D-α-tocopheryl polyethylene glycol 1000 succinate (TPGS)"}


def sentence_case(t):
    """Crossref titles are sometimes in title case; convert to sentence case unless a word looks like a name, acronym or identifier."""
    ws = t.split(" ")
    caps = sum(1 for w in ws[1:] if len(w) > 3 and w[0].isupper() and w[1:].islower())
    if caps < 0.5 * max(1, sum(1 for w in ws[1:] if len(w) > 3)):
        return t
    out = [ws[0]]; prev = ws[0]
    for w in ws[1:]:
        if w[0].isupper() and w[1:].islower() and not prev.endswith(":"):
            w = w.lower()
        out.append(w); prev = w
    r = " ".join(out)
    for a, b in (("p-glycoprotein", "P-glycoprotein"), ("bayesian", "Bayesian"), ("Inhibition of P-Glycoprotein by D-α-Tocopheryl Polyethylene Glycol 1000 Succinate", "Inhibition of P-glycoprotein by D-α-tocopheryl polyethylene glycol 1000 succinate")):
        r = r.replace(a, b)
    return r


def build():
    lib = {}
    for k, m in META.items():
        au = m["authors"]; yr = YEAR.get(k, m["year"])
        fam = [a.split(",")[0] for a in au]
        if len(au) == 1: short = fam[0]
        elif len(au) == 2: short = "%s and %s" % (fam[0], fam[1])
        else: short = "%s et al." % fam[0]
        alist = ", ".join(au[:6]) + (", et al." if len(au) > 6 else "")
        jr = ABBR.get(m["journal"], m["journal"]); m = dict(m, title=TITLE_OVERRIDE.get(k, sentence_case(m["title"])))
        vol = m["volume"]; iss = ("(%s)" % m["issue"]) if m["issue"] and m["journal"] not in ("Pharmaceutics", "Scientific Reports") else ""
        pg = m["pages"].replace("-", "–")
        full = "%s, %d. %s. *%s* %s%s, %s. https://doi.org/%s" % (alist, yr, m["title"], jr, vol, iss, pg, m["doi"]) if pg else "%s, %d. %s. *%s* %s. https://doi.org/%s" % (alist, yr, m["title"], jr, vol, m["doi"])
        lib[k] = (short, yr, fam[0] + " %d" % yr, full)
    return lib

LIB = build()
if __name__ == "__main__":
    for k, v in LIB.items(): print(k, "|", v[0], v[1]); print("   ", v[3])
