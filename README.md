# EPSR_hackathon2026

A Python port of EPSR, subject of the 2026 "Reinvent The Wheel" hackathon run by Science Works.

* [python/](python/README.md): the Python package `epsr` (install, run, examples on real experimental data).
* [rust/](rust/README.md): the Rust translation of the Fortran EPSR engine that the package drives.
* [validation/](validation/README.md): water, silica and Cu-BTC refined against measured neutron data by the Fortran, the Rust CLI and the Python wrapper.
* [mellor/](mellor/docs/integration.md): the LLM agent that reads measured data, proposes a model system and drives a refinement.
* [webapp/](webapp/README.md): the web application around Mellor (uploads, accounts, a job queue, seaborn figures, PDF reports, and publishing figures to an Overleaf project).
* `EPSR-test/`, `test_problems/`: the test case and verification cases; `experimental_data/`: the Dissolve examples; `fortran_src/`: the original Fortran.

