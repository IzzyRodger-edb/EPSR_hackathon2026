# EPSR-test-ex32 — `mixato` (file selection, error paths)

**Program:** `mixato`.

**Inputs:** `DMSO.ato`, `Water.ato` (= `EPSR-test/WaterMolecule.ato` renamed, so its
template `WaterMolecule.ato` does not exist), a Linux `system_commands.txt`
(`system_ls ls`, `system_termination /`), `mixato.stdin`.

**Command:** `mixato 2 < mixato.stdin` (number of files on the command line).

**Exercises** `search_for_files` with a real directory listing (`ls ./*.ato > dirlist.txt`):
`u` (go back) on the first entry (re-positioning at the end of the list), `y` accept,
`n` rejections with wrap-around (rewind), `e` (exit list) and typing a name; a
non-existent file (`NoSuchFile` → "Problems with .ato file ... Try again");
the missing-template branch (template renamed to the file name `Water`);
output name with an extension (`mix2.out` → `mix2.ato`).

**Expected outputs:** `mix2.ato`, `dirlist.txt` (directory listing; its order is that of
`ls` in the current locale — the answers are chosen so that the result does not depend on it).

**Verification:** ALL IDENTICAL.
