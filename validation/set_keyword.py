#!/usr/bin/env python3
"""Set the value of a keyword in an EPSR keyword file (.EPSR.inp, .pcof, .NWTS.dat, ...).

usage: set_keyword.py FILE KEYWORD VALUE [VALUE ...]

This is the edit a user makes in the EPSR shell (or a text editor) between runs,
e.g. switching on the empirical potential with `set_keyword.py run.pcof ereq 3.0`.
Only the keyword block (up to the first line that is just `q`) is searched.  The
line is rewritten as `<keyword padded to 12><value>               <description>`,
the layout EPSR itself writes, and the program re-reads it with list-directed input.
"""
import sys


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        return 2
    path, kw, value = sys.argv[1], sys.argv[2], " ".join(sys.argv[3:])
    with open(path, encoding="latin-1") as f:
        lines = f.read().split("\n")
    for i, line in enumerate(lines):
        if line.strip() == "q":
            break
        if line[:12].strip() == kw:
            desc = line.rstrip().split(" " * 15)[-1] if " " * 15 in line.rstrip() else ""
            lines[i] = f"{kw:<12}{value}               {desc}"
            with open(path, "w", encoding="latin-1") as f:
                f.write("\n".join(lines))
            return 0
    print(f"set_keyword.py: keyword {kw!r} not found in {path}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
