# Presentation assets

Movies and still images of the two EPSR examples, for slides. Rebuild them
with:

```bash
cd visualisation/python_visualisation_toolkit
uv sync --extra viewer
uv run python presentation_assets/make_presentation_assets.py            # all (~10 min)
uv run python presentation_assets/make_presentation_assets.py --stills   # stills only (~1 min)
```

| file | content |
|---|---|
| `H2O_config_movie.mp4` | water: 3D configuration (10 Å slab) + partial g(r), iterations 1 → 300 (one frame per 2 iterations, 12 fps, 1920×1080, H.264). The g(r) panel compares with the **experimental partials**: experiment solid, EPSR the same colour but translucent |
| `H2O_config_movie_2.mp4` | as `H2O_config_movie.mp4`, but the experimental partials fade out over iterations 1 → 100 (EPSR curves become opaque), the final distribution is held for 3 s, then the experiment fades back in. The R² of each EPSR partial vs experiment (0–10 Å, 4 d.p.) is shown under the legend on every frame. Made by `make_H2O_movie_2.py` (~1 min) |
| `SiO2_config_movie.mp4` | silica: 3D configuration (12 Å slab, Si coloured by number of O neighbours) + calculated partial g(r), same settings |
| `H2O_config_iter100.png`, `H2O_config_iter300.png` | water movie frames at the end of the reference-potential stage and at the end |
| `SiO2_config_iter100.png`, `SiO2_config_iter300.png` | the same for silica |
| `H2O_partials_vs_experiment.png` | g_OO, g_OH, g_HH at iteration 300 vs experiment, one panel per pair (matplotlib, 200 dpi) |
| `fluid_configuration.xyz` | example fluid: 400 atoms of a 2D Lennard-Jones liquid (ρ* = 0.70, T* = 1.0, reduced units, σ = 1) from a Metropolis Monte Carlo run; extended XYZ, z = 0 |
| `fluid_configuration_gr.dat` | its simulated g(r) (peaks at 1.14, 2.21, 3.29 σ) |
| `fluid_configuration.png` | plain rendering of that configuration |
| `fluid_gr_cartoon.png` | cartoon g(r) of a typical fluid, annotated: excluded core, 1st and 2nd shells, g(r) → 1. Peak positions taken from the simulation |
| `fluid_config_cartoon.png` | cartoon of the configuration seen from one atom: the 1st-shell neighbours, the 2nd shell, and the r → r + dr shell that g(r) counts |
| `fluid_jiggle.gif` | the example fluid's atoms jiggling: the Monte Carlo continued from the saved configuration, one sweep per frame (120 frames, ~8 s loop, 640×640) |
| `make_fluid_cartoons.py` | the script that makes the six fluid files (`--rerun` repeats the ~1 min simulation; otherwise it reuses the saved configuration) |
| `experimental_data/` | the experimental files and their source (`SOURCE.md`) |
| `make_presentation_assets.py` | the script that makes all of the above |

The MP4s are full-HD versions of the examples' `config_movie.gif`. They are
rendered directly from the trajectories rather than converted from the GIFs, so
they have full colour depth instead of the GIF's 256-colour palette.

## Experimental partial pair distribution functions

The water partials are from A. K. Soper, *Chem. Phys.* **258** (2000) 121,
downloaded from the [ISIS Disordered Materials Database](http://wwwisis2.isis.rl.ac.uk/Disordered/Database/DBMain.htm).
They are the 298 K, ambient-pressure set ("case (a)"; see
`experimental_data/SOURCE.md` for how the case was identified). The simulation
is at 300 K and 0.1 atoms/Å³.

- The experimental partials include intramolecular correlations, so the EPSR
  curves are inter + intramolecular g(r).
- The y-axis is scaled to the intermolecular structure, so the intramolecular
  O–H bond peak at ~1 Å goes off the top of the panel.
- **Silica has no experimental partials in the database.** It provides only the
  total correlation function T(r) of vitreous SiO₂, so the silica movie shows
  the calculated partials alone.

What the comparison shows at iteration 300 (peak positions and heights, EPSR
vs experiment):

| feature | EPSR | experiment |
|---|---|---|
| O–O first peak | 2.79 Å, 2.5 | 2.73 Å, 2.75 |
| O–O first minimum | 3.45 Å, 0.95 | 3.36 Å, 0.78 |
| O–O second shell (tetrahedral) | no distinct peak near 4.5 Å | 4.5 Å, 1.16 |
| O···H hydrogen-bond peak | 1.83 Å, 0.95 | 1.77 Å, 1.12 |
| H–H first intermolecular peak | 2.43 Å, 1.18 | 2.34 Å, 1.28 |

- **Peak positions agree** to within ~0.1 Å, and the H–H partial agrees well.
- **The EPSR first shells are lower and broader** than the experimental ones
  (O–O and the O···H hydrogen bond).
- **The tetrahedral second O–O shell at 4.5 Å is essentially missing**, so the
  model water is less ordered than the reference. Likely contributors: the short
  refinement (300 iterations) with a small EP amplitude (`ereq` = 3 kJ/mol); and
  different input data (the SANDALS data sets of the Dissolve example, not
  those of the 2000 paper). This was not tested further.
- The EPSR values above were read from a 5-point running mean, to smooth the
  noise.

Note that these "experimental" partials are not a direct measurement. Soper
extracted them from H₂O/D₂O/HDO diffraction data with an EPSR-type
refinement, so they are the standard reference for water but are themselves
model-assisted.

## In the toolkit

These additions are available to any script:

- `epsrvis.read_isis_gr(path)` / `epsrvis.load_reference_partials({pair: file})`:
  read ISIS-database partial g(r) files.
- `ConfigurationViewer(..., reference=..., calc_opacity=0.4)`: overlay reference
  partials in the g(r) panel, with the calculated g(r) translucent.
  `window_size=(1920, 1080)` gives full-HD frames; text scales with the height.
- `ConfigurationViewer.movie("x.mp4", quality=9)`: H.264 MP4 (needs
  `imageio-ffmpeg`, part of the `viewer` extra).
- `plots.partials_vs_reference(traj, reference)`: the static comparison figure.
