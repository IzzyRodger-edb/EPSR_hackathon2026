#!/usr/bin/env python3
"""Render the one-page summary PDF of every verification case.

usage: make_case_pdfs.py [--root REPO_ROOT] [--results REGRESSION_LOG] [N ...]

For each `test_problems/<system>/EPSR-test-ex{N}` folder of the repository, the summary
`rust/docs/case_summaries/ex{N}.json` is rendered to
`EPSR-test-ex{N}/EPSR-test-ex{N}.pdf` (A4, one page).  The number of chained
passes and the programs run are taken from the folder's `commands.txt`; the
verification result is taken from a log of `run_all_ex_tests.sh` if given
(lines `EPSR-test-exN: ALL IDENTICAL (k pass(es))`).

Summary JSON fields (strings may use `backticks` for code):
  title        short title of the case
  purpose      what feature / code path the case exercises (1-3 sentences)
  fortran      list of Fortran files/routines exercised
  rust         list of Rust files/functions of the port
  system       the physical system simulated / processed
  inputs       how the inputs were made
  verification how the accuracy was checked (beyond the standard comparison)
  notes        optional: Fortran quirks reproduced, limitations
  refs         list of keys of case_summaries/references.json; cite them in the
               text as [@key] or [@key1,@key2] (numbered in order of first use)
  extra_refs   optional list of additional full references

No third-party packages: the PDF is written directly (standard Type 1 fonts
Helvetica / Helvetica-Bold / Courier, WinAnsiEncoding), with text wrapped
using the Adobe font metrics and the font size reduced until the page fits.
"""
import json
import os
import re
import sys
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- font metrics
_HELV = dict(zip(range(32, 127), [
    278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278,
    556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584, 556,
    1015, 667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722, 778,
    667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469, 556,
    333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556,
    556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584]))
_HELVB = dict(zip(range(32, 127), [
    278, 333, 474, 556, 556, 889, 722, 238, 333, 333, 389, 584, 278, 333, 278, 278,
    556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 333, 333, 584, 584, 584, 611,
    975, 722, 722, 722, 722, 667, 611, 778, 722, 278, 556, 722, 611, 833, 722, 778,
    667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 333, 278, 333, 584, 556,
    333, 556, 611, 556, 611, 556, 333, 611, 611, 278, 278, 556, 278, 889, 611, 611,
    611, 611, 389, 556, 333, 611, 556, 778, 556, 556, 500, 389, 280, 389, 584]))
# cp1252 characters outside ASCII that the summaries may use
_EXTRA = {0xC5: (667, 722), 0xB0: (400, 400), 0xD7: (584, 584), 0x96: (556, 556), 0x97: (1000, 1000),
          0xB1: (584, 584), 0xB2: (333, 333), 0xB3: (333, 333), 0xB7: (278, 278), 0xE9: (556, 556),
          0xF6: (556, 611), 0xFC: (556, 611), 0xE4: (556, 556), 0xB5: (556, 611), 0x92: (222, 278),
          0x91: (222, 278), 0x93: (333, 500), 0x94: (333, 500), 0x85: (1000, 1000), 0xD6: (778, 778),
          0xBD: (834, 834), 0x95: (350, 350), 0xC3: (667, 722), 0xE5: (556, 556)}
_SUBST = {"ρ": "rho", "θ": "theta", "ω": "omega", "φ": "phi", "ψ": "psi", "Δ": "Delta", "α": "alpha",
          "β": "beta", "γ": "gamma", "λ": "lambda", "σ": "sigma", "ε": "epsilon", "π": "pi", "χ": "chi",
          "μ": "µ", "≈": "~", "≥": ">=", "≤": "<=", "→": "->", "←": "<-", "⇒": "=>", "∝": "~",
          "⟨": "<", "⟩": ">", "∞": "inf", "√": "sqrt", "≠": "/=", "·": "·", "−": "-", "…": "...",
          "Å": "Å", "‘": "'", "’": "'", "“": '"', "”": '"', "×": "×"}


def clean(s):
    return "".join(_SUBST.get(c, c) for c in s)


def enc(s):
    return s.encode("cp1252", errors="replace")


def width(s, font, size):
    w = 0
    for b in enc(s):
        if font == "C":
            w += 600
        elif b in _HELV:
            w += (_HELVB if font == "B" else _HELV)[b]
        elif b in _EXTRA:
            w += _EXTRA[b][1 if font == "B" else 0]
        else:
            w += 556
    return w * size / 1000.0


# ---------------------------------------------------------------- layout
def runs(text, font="R"):
    """Split text into (word, font, space_before) tokens; `code` -> Courier."""
    out = []
    parts = re.split(r"(`[^`]*`)", clean(text))
    pending_space = False
    for p in parts:
        if not p:
            continue
        if p.startswith("`") and p.endswith("`"):
            words = p[1:-1].split(" ")
            f = "C"
        else:
            words = re.split(r"(\s+)", p)
            f = font
        first = True
        for w in words:
            if f != "C" and w and w.isspace():
                pending_space = True
                continue
            if f == "C" and not first:
                pending_space = True
            first = False
            if w == "":
                continue
            out.append((w, f, pending_space))
            pending_space = False
        if f != "C" and p and p[-1].isspace():
            pending_space = True
    return out


class Page:
    W, H = 595.28, 841.89
    ML, MR, MT, MB = 48, 48, 46, 40

    def __init__(self, size):
        self.size = size
        self.ops = []
        self.y = self.H - self.MT
        self.overflow = False

    def need(self, h):
        if self.y - h < self.MB:
            self.overflow = True
        self.y -= h

    def text(self, x, y, s, font, size, color=None):
        name = {"R": "F1", "B": "F2", "C": "F3", "I": "F4"}[font]
        col = "%.3f %.3f %.3f rg " % color if color else ""
        esc = enc(s).replace(b"\\", b"\\\\").replace(b"(", b"\\(").replace(b")", b"\\)")
        self.ops.append(b"BT " + col.encode() + ("/%s %.2f Tf %.2f %.2f Td (" % (name, size, x, y)).encode()
                        + esc + b") Tj ET")

    def rule(self, y, gray=0.6, w=0.6):
        self.ops.append(("%.2f G %.2f w %.2f %.2f m %.2f %.2f l S 0 G" % (gray, w, self.ML, y, self.W - self.MR, y)).encode())

    def band(self, y, h, rgb):
        self.ops.append(("%.3f %.3f %.3f rg %.2f %.2f %.2f %.2f re f 0 g" % (rgb + (self.ML - 6, y, self.W - self.ML - self.MR + 12, h))).encode())

    def para(self, text, size=None, indent=0, bullet=None, font="R", lead=1.22, color=None):
        size = size or self.size
        x0 = self.ML + indent
        maxw = self.W - self.MR - x0
        toks = runs(text, font)
        # tokens not preceded by a space (e.g. punctuation after a code span)
        # are kept on the same line as the token before them
        units = []
        for w, f, sp in toks:
            if units and not sp:
                units[-1][1].append((w, f))
            else:
                units.append((sp, [(w, f)]))
        line, lw = [], 0.0
        lines = []
        for sp, parts in units:
            items = [(w, f, size * (0.92 if f == "C" else 1.0)) for w, f in parts]
            uw = sum(width(w, f, fs) for w, f, fs in items)
            spw = width(" ", "R", size) if (sp and line) else 0.0
            if line and lw + spw + uw > maxw:
                lines.append(line)
                line, lw, spw = [], 0.0, 0.0
            for k, (w, f, fs) in enumerate(items):
                line.append((w, f, fs, spw if k == 0 else 0.0))
            lw += spw + uw
        if line:
            lines.append(line)
        for i, ln in enumerate(lines):
            self.need(size * lead)
            if i == 0 and bullet:
                self.text(x0 - (17 if len(bullet) > 1 else 9), self.y, bullet, "R", size, color)
            x = x0
            for w, f, fs, spw in ln:
                x += spw
                self.text(x, self.y, w, f, fs, color)
                x += width(w, f, fs)

    def heading(self, s):
        self.need(self.size * 0.55)
        self.need(self.size * 1.35)
        self.text(self.ML, self.y, s, "B", self.size * 1.08, (0.10, 0.25, 0.50))
        self.rule(self.y - 2.5, 0.75, 0.4)
        self.y -= 2.0


def number_citations(summ, bib):
    """Replace `[@key]` / `[@a,@b]` citations by numbers, in order of first use;
    keys listed in `refs` but not cited follow."""
    order = []
    fields = ["purpose", "system", "inputs", "verification", "notes"]
    pat = re.compile(r"\[(@[A-Za-z0-9_]+(?:\s*,\s*@[A-Za-z0-9_]+)*)\]")
    for f in fields:
        for m in pat.finditer(summ.get(f, "") or ""):
            for k in re.findall(r"@([A-Za-z0-9_]+)", m.group(1)):
                if k not in order:
                    order.append(k)
    for k in summ.get("refs", []):
        if k not in order:
            order.append(k)
    bad = [k for k in order if k not in bib]
    num = {k: i + 1 for i, k in enumerate(order)}
    out = dict(summ)
    for f in fields:
        if out.get(f):
            out[f] = pat.sub(lambda m: "[" + ",".join(str(num.get(k, "?")) for k in re.findall(r"@([A-Za-z0-9_]+)", m.group(1))) + "]", out[f])
    out["refs"] = order
    return out, bad


def build(case, summ, meta, bib, size):
    pg = Page(size)
    # title band
    pg.band(pg.H - pg.MT - 26, 34, (0.90, 0.93, 0.98))
    pg.y -= 2
    pg.text(pg.ML, pg.y - 6, f"EPSR-test-ex{case}", "B", 16, (0.10, 0.25, 0.50))
    tw = width(f"EPSR-test-ex{case}  ", "B", 16)
    title = clean(summ["title"]).replace("`", "")
    tsize = 12.5
    while width(title, "B", tsize) > pg.W - pg.ML - pg.MR - tw and tsize > 8:
        tsize -= 0.25
    pg.text(pg.ML + tw, pg.y - 6, title, "B", tsize)
    pg.y -= 30
    info = (f"Program(s): {meta['programs']}   |   chained passes: {meta['passes']}   |   "
            f"verification: {meta['result']}")
    pg.para(info, size=size * 0.92, color=(0.25, 0.25, 0.25))
    pg.heading("What this case tests")
    pg.para(summ["purpose"])
    pg.heading("Code exercised")
    pg.para("Fortran (original, fortran_src/):", font="B")
    for s in summ["fortran"]:
        pg.para(s, indent=14, bullet="•" if False else "-")
    pg.para("Rust (port, rust/epsr-core/src/):", font="B")
    for s in summ["rust"]:
        pg.para(s, indent=14, bullet="-")
    pg.heading("Physical system")
    pg.para(summ["system"])
    pg.heading("Inputs")
    pg.para(summ["inputs"])
    pg.heading("How the accuracy was verified")
    pg.para(meta["standard"])
    if summ.get("verification"):
        pg.para(summ["verification"])
    if summ.get("notes"):
        pg.heading("Notes")
        pg.para(summ["notes"])
    refs = [bib[k] for k in summ.get("refs", [])] + list(summ.get("extra_refs", []))
    if refs:
        pg.heading("References")
        for i, r in enumerate(refs, 1):
            pg.para(r, size=size * 0.85, indent=20, bullet=f"[{i}]", lead=1.15)
    # footer
    pg.rule(pg.MB - 8, 0.75, 0.4)
    pg.text(pg.ML, pg.MB - 19, f"EPSR Fortran-to-Rust port - verification case ex{case}. "
            f"Source: rust/docs/case_summaries/ex{case}.json, rendered by rust/tools/make_case_pdfs.py.",
            "R", 6.8, (0.4, 0.4, 0.4))
    return pg


def pdf_bytes(pg, title):
    content = zlib.compress(b"\n".join(pg.ops))
    objs = []
    objs.append(b"<< /Type /Catalog /Pages 2 0 R >>")
    objs.append(b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>")
    objs.append(("<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.2f %.2f] /Resources << /Font << "
                 "/F1 5 0 R /F2 6 0 R /F3 7 0 R /F4 8 0 R >> >> /Contents 4 0 R >>" % (pg.W, pg.H)).encode())
    objs.append(b"<< /Length %d /Filter /FlateDecode >>\nstream\n" % len(content) + content + b"\nendstream")
    for base in ("Helvetica", "Helvetica-Bold", "Courier", "Helvetica-Oblique"):
        objs.append(("<< /Type /Font /Subtype /Type1 /BaseFont /%s /Encoding /WinAnsiEncoding >>" % base).encode())
    t = enc(title).replace(b"\\", b"\\\\").replace(b"(", b"\\(").replace(b")", b"\\)")
    objs.append(b"<< /Title (" + t + b") /Producer (make_case_pdfs.py) >>")
    out = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offs = []
    for i, o in enumerate(objs, 1):
        offs.append(len(out))
        out += b"%d 0 obj\n" % i + o + b"\nendobj\n"
    xref = len(out)
    out += b"xref\n0 %d\n0000000000 65535 f \n" % (len(objs) + 1)
    for o in offs:
        out += b"%010d 00000 n \n" % o
    out += b"trailer\n<< /Size %d /Root 1 0 R /Info %d 0 R >>\nstartxref\n%d\n%%%%EOF\n" % (len(objs) + 1, len(objs), xref)
    return bytes(out)


def case_meta(d, n, results):
    progs, rep, cmds = [], 1, 0
    for line in open(os.path.join(d, "commands.txt")):
        w = line.split("#")[0].split()
        if not w:
            continue
        if w[0] == "repeat":
            rep = int(w[1])
        else:
            cmds += 1
            if w[0] not in progs:
                progs.append(w[0])
    res = results.get(n)
    if res is None:
        result = "see README.md"
        tail = ""
    else:
        ok, k = res
        np_ = "1 pass" if k == 1 else f"all {k} passes"
        result = (f"ALL IDENTICAL ({np_})" if ok else "DIFFERENCES - see README.md")
        tail = (f" Result of the final regression run of the complete port: {'every file identical in ' + np_ if ok else 'differences found'}.")
    patch = ""
    if {"readxyz", "readpdb"} & set(progs):
        patch = (" `readxyz`/`readpdb` crash (SIGSEGV) as compiled from the original source, so the reference build adds an "
                 "explicit INTERFACE block for `ctype` (`rust/tools/build_patched_readxyz_readpdb.sh`); the port behaves "
                 "like the patched program.")
    elif {"sharm", "sdf"} & set(progs) or n == 122:
        patch = (" A run-time FORMAT in the original source (`sharm_rev2.f`/`sdf_rev2.f`, and `sharm_diffuse` for `epsr`) lacks "
                 "a comma after the Hollerith `1h(` and is rejected by current gfortran, so the reference build inserts that "
                 "comma (originals kept as `*.orig`); the port behaves like the patched program.")
    masks = ("only the date line of `.out` and the wall-clock columns of `fort.35` are masked" if "epsr" in progs
             else "no masking except where noted below")
    calls = "1 program call" if cmds == 1 else f"{cmds} program calls"
    chain = ("run once" if rep == 1 else f"repeated {rep} times, each pass reading the files written by the previous one")
    std = (f"`commands.txt` ({calls}, {chain}) was run with the Fortran reference build (`rust/tools/build_fortran_reference.sh`: the original "
           f"sources, gfortran -O2 -fno-tree-vectorize, single thread) and with the Rust `epsr-rs` on identical copies of "
           f"the inputs, and every file written was compared byte for byte (`rust/tools/run_ex_test.sh`, "
           f"`compare_outputs.py`; {masks}). The Fortran output is stored in `fortran_reference.tar.gz`." + patch + tail)
    return {"programs": ", ".join(progs), "passes": rep, "result": result, "standard": std}


def main():
    args = sys.argv[1:]
    root = os.path.normpath(os.path.join(HERE, "..", ".."))
    results_file = None
    sel = []
    i = 0
    while i < len(args):
        if args[i] == "--root":
            root = args[i + 1]; i += 2
        elif args[i] == "--results":
            results_file = args[i + 1]; i += 2
        else:
            sel.append(int(args[i])); i += 1
    sdir = os.path.join(HERE, "..", "docs", "case_summaries")
    bib = json.load(open(os.path.join(sdir, "references.json"), encoding="utf-8"))
    results = {}
    if results_file:
        for line in open(results_file):
            m = re.match(r"EPSR-test-ex(\d+): (ALL IDENTICAL|DIFFERENCES)\D*(\d+)", line)
            if m:
                k = int(m.group(3)) if m.group(2) == "ALL IDENTICAL" else 0
                results[int(m.group(1))] = (m.group(2) == "ALL IDENTICAL", k)
    tp = os.path.join(root, "test_problems")
    folders = {}
    for system in sorted(os.listdir(tp)):
        if not os.path.isdir(os.path.join(tp, system)):
            continue
        for x in os.listdir(os.path.join(tp, system)):
            m = re.fullmatch(r"EPSR-test-ex(\d+)", x)
            if m:
                folders[int(m.group(1))] = os.path.join(tp, system, x)
    cases = sorted(folders)
    if sel:
        cases = [c for c in cases if c in sel]
    nfail = 0
    for n in cases:
        d = folders[n]
        jf = os.path.join(sdir, f"ex{n}.json")
        if not os.path.exists(jf):
            print(f"ex{n}: no summary {jf}"); nfail += 1
            continue
        summ = json.load(open(jf, encoding="utf-8"))
        summ, missing = number_citations(summ, bib)
        if missing:
            print(f"ex{n}: unknown reference keys {missing}"); nfail += 1
            continue
        meta = case_meta(d, n, results)
        size = 9.6
        while True:
            pg = build(n, summ, meta, bib, size)
            if not pg.overflow or size < 6.6:
                break
            size -= 0.15
        if pg.overflow:
            print(f"ex{n}: WARNING text does not fit on one page"); nfail += 1
        with open(os.path.join(d, f"EPSR-test-ex{n}.pdf"), "wb") as f:
            f.write(pdf_bytes(pg, f"EPSR-test-ex{n}: {clean(summ['title']).replace('`', '')}"))
        print(f"ex{n}: EPSR-test-ex{n}.pdf (font {size:.2f} pt)")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
