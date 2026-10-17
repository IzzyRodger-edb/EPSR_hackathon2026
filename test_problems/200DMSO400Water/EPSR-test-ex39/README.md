# EPSR-test-ex39 — box-building chain: mixato → makelatticeato → randomise → fmole

**Programs:** `mixato`, `makelatticeato`, `randomise`, `fmole`, with a Linux
`system_commands.txt` (`system_ls ls`, `system_termination /`) so that every file prompt
goes through the real directory listing (answered with `e` + typed name, which does not
depend on the `ls` sort order).

**Commands:**
1. `mixato < mixato.stdin` — 4 DMSO + 8 water, density 0.1 → `mixbox.ato`.
2. `makelatticeato mixbox 2 2 2 big n` — 96-molecule cubic box.
3. `randomise big`.
4. `fmole < fmole.stdin` — file chosen via the listing prompt, 10 shakes.

**Expected outputs:** `mixbox.ato`, `big.ato`, `dirlist.txt` (listing of `./*.ato`).

**Verification:** ALL IDENTICAL.
