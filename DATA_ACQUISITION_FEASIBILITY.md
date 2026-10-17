# Feasibility: an agent that acquires experimental data for EPSR

**Status:** feasibility assessment and specification. No code has been written. The
acquisition agent is to be built separately; this document is what it should be built
from.

**Question asked:** can an agent fetch experimental neutron-scattering data from
facilities such as ISIS (STFC, UK) or Oak Ridge National Laboratory (SNS/HFIR, US) and
feed it to this EPSR port?

**Scope fixed by the requester:** pluggable sources (one driver per repository, public
and credentialed); credentials optional and modular, so the agent works without any;
**already-reduced data only** — driving Gudrun or Mantid to reduce raw data is out of
scope.

Every external claim below is marked **[CONFIRMED]** (verified against official
documentation or by exercising the endpoint), **[PARTIAL]** (official-adjacent but not
fully verified) or **[UNCONFIRMED]** (open question for the implementer). Internal
claims cite `file:line` in this repo.

---

## 1. Verdict

Yes, but **not in the shape the question implies**, and the reason matters more than the
answer.

| Source | Reduced, EPSR-ready data? | Credentials | Verdict |
|---|---|---|---|
| **Public DOI'd archives** (Bath, Cambridge Apollo, Dissolve examples) | **Yes — actual `.mint01`** | None | **Primary source. Start here.** |
| **ISIS catalogue** (DataGateway/ICAT API) | **No — RAW/NeXus only** | STFC User Office account | API works, but cannot deliver in-scope data |
| **ORNL** (ONCat API) | No — event NeXus; no `.mint01` path exists | Anonymous for metadata; issued credentials for bytes | Discovery useful; conversion unsolved |
| **ILL** (D4c) | Unknown holdings | Registration | Better second facility than ORNL (Gudrun supports D4c) |

Three findings drive that table.

**Finding 1 — the format is trivial; this is not the hard part.** EPSR's `nrtype 5`
reader wants plain text: header lines whose first whitespace token contains `#`, then
three columns `Q S(Q) E(Q)`. Writing a conforming file is an afternoon's work (§2.1).

**Finding 2 — the obvious plan does not work: ISIS does not publish reduced data.** The
ISIS Data Policy (May 2025) distinguishes *raw data* (s.2.1), *facility-generated reduced
data* (s.2.2, retained a minimum of six months only) and *result data* (s.2.4 — the
output of manual reduction, which is what `.mint01` is). Long-term curation and public
DOI release apply to **raw RAW/NeXus**; the policy states long-term storage of result
data via eData is for STFC staff and that "ISIS currently does not provide any such
facility for non-STFC staff" (s.4.2.3). **[CONFIRMED]**
→ So "query the ISIS portal, download, run EPSR" is not a workflow that exists. You would
be downloading raw detector data and reducing it yourself with Gudrun — which the
requester placed out of scope. The ISIS catalogue API is real and reachable, but it is a
*raw-data* API.

**Finding 3 — the real bottleneck is sample metadata, not bytes.** A scattering file
cannot start a refinement. EPSR also needs the composition, number density, isotopic
substitution scheme and a reference force field — none of which exists in a data file or
in any facility catalogue (§2.2). For the three systems already in this repo, that
metadata came from the Dissolve `input.txt` files, not from the data.

**Consequence for the design:** the highest-value driver is a *public-archive* driver,
because a few archived records contain reduced data **and** the metadata — one of them
(§3.1) is a complete EPSR project including `.mol`, `.ato`, `.NWTS` and `.EPSR.inp`.
Facility catalogue APIs are worth implementing for **discovery and provenance**, and for
the day raw→reduced comes into scope.

---

## 2. What EPSR actually needs

### 2.1 The data-file contract (easy)

Reader: [`rust/epsr-core/src/datafile.rs:81-171`](rust/epsr-core/src/datafile.rs#L81-L171).
Accepted types, from [`fortran_src/epsr_standalone_rev3.f:334-340`](fortran_src/epsr_standalone_rev3.f#L334-L340)
and the Rust module doc:

| `nrtype` | Format |
|---|---|
| 1 | `Q S(Q)`, list-directed, one pair per line, errors set to 0 |
| 2 | `Q S(Q) E(Q)`, list-directed |
| 3, 4 | GENIE "SHOW DATA": `N Q S(Q) E(Q)` as `format(i7,3g14.7)` after 4 header lines (3 = histogram, 4 = point mode) |
| **5** | **Gudrun `.mint01`: `Q S(Q) E(Q)`, list-directed, lines whose first token contains `#` are skipped (histogram)** |

For `nrtype 5`, the rules an emitter must satisfy:

- **Headers**: any number of leading lines, each with `#` **in its first whitespace
  token** (`# title` works; a line beginning `  x # y` does not skip). No fixed count.
- **Columns**: three whitespace-separated floats. Any float syntax — read list-directed,
  so exact column positions are irrelevant despite the real files being fixed-width.
- **Q is a histogram bin edge.** EPSR replaces each Q with the midpoint to the next and
  **drops the last point** ([`epsr.rs:502-509`](rust/epsr-core/src/epsr.rs#L502-L509)).
  Emit a final throwaway row; all four real files in this repo end with
  `Q 0.0000000E+00 0.0000000E+00`.
- **Limits**: ≤ 5000 points (`MQD`, [`state.rs:31`](rust/epsr-core/src/state.rs#L31)).
- **Placement**: the path is built as `directory + name`
  ([`epsr.rs:198-220`](rust/epsr-core/src/epsr.rs#L198-L220)) and truncated at the first
  blank, so the file must sit **directly in the run folder** with a **blank-free name**.
  Dissolve-style `data/` subdirectories do not work.
- **Useful Q range**: data is interpolated onto EPSR's internal grid `q[i] = qstep*(i-1)`,
  `i = 1..nq`. With the validation settings (`nq 600`, `qstep 0.05`) anything above
  ≈30 Å⁻¹ **is read and silently discarded**. Worth warning about in the validator.
- The extension is **never checked**. `.mint01` is convention. The name before the first
  `.` (max 28 chars) becomes the column label in the `.u01`/`.v01` outputs
  ([`epsr.rs:1022-1024`](rust/epsr-core/src/epsr.rs#L1022-L1024)).

A working 8-line reader to reuse as the reference implementation already exists:
[`test_problems/200DMSO400Water/EPSR-test-ex95/scripts/make_data.py`](test_problems/200DMSO400Water/EPSR-test-ex95/scripts/make_data.py),
which converts `.mint01` to the `nrtype` 1/3/4 forms.

### 2.2 The metadata gap (hard — this is the real work)

To refine against a new dataset, EPSR additionally needs:

| Need | File | Derivable from scattering data? |
|---|---|---|
| Molecular geometry + **reference force field** (LJ ε/σ, mass, charge, symbol) | `.mol` per species, e.g. [`validation/water/Water.mol`](validation/water/Water.mol) | **No** |
| Species populations + **number density** (atoms/Å³) | `mixato` stdin answers | **No** |
| **Isotopic substitution** per sample (`abundances`, `iexchange`), matching `.ato` atom labels exactly | `.NWTS.dat` per sample → run `nwts` → `.NWTStot.wts` | **No** |
| Fit settings and dataset wiring | `.EPSR.inp` — EPSR never creates one | Partly (conventions) |

The `.EPSR.inp` `data n` block is **positional**: the eight keywords are read in fixed
order — `datafile`, `wtsfile`, `nrtype`, `rshmin`, `szeros`, `tweak`, `efilereq`,
`porodadd` — with the keyword in columns 1–12, 15 spaces before the description, and the
block closed by a line that is just `q`. `ndata` must be incremented to match.
([`epsr.rs:198-220`](rust/epsr-core/src/epsr.rs#L198-L220); layout per
[`validation/set_keyword.py`](validation/set_keyword.py).)

The full setup chain (`makemole` → `mixato` → `randomise` → `fmole` → `nwts` → write
`.EPSR.inp` → first scratch run → equilibrate → enable EP → refine) is already documented
in [`validation/PYTHON_WRAPPER_WORKFLOW.md`](validation/PYTHON_WRAPPER_WORKFLOW.md) §1 and
exercised by [`validation/*/steps.txt`](validation/water/steps.txt). **The agent should
target that existing chain, not reinvent it.**

**Design implication:** keep the metadata resolver a *separate, explicit, human-reviewable*
stage. Guessing a force field from a diffraction pattern is not possible; presenting a
candidate for a human to confirm is. Where an archive ships the metadata (§3.1), ingest
it rather than reconstruct it.

---

## 3. Source-by-source assessment

### 3.1 Public DOI'd archives — the primary source **[CONFIRMED]**

These need no credentials, carry explicit licences, and contain real reduced data.

**Bath Research Data Archive, DOI [`10.15125/BATH-00359`](https://researchdata.bath.ac.uk/359/)** — CC-BY-4.0, ~654 kB.
Reduced, background-corrected wide-Q neutron data **in `.mint01` format**, processed with
GudrunN, hydrogen inelastic scattering subtracted, described by the depositors as
"intended as input for EPSR modelling". Instruments NIMROD and SANDALS; deep eutectic
solvent + water mixtures across multiple H/D labellings and dilutions. Contributor
Daniel Bowron (ISIS). A sibling record exists at
[researchdata.bath.ac.uk/401](https://researchdata.bath.ac.uk/401/) **[PARTIAL — not individually verified]**.

**Cambridge Apollo, DOI [`10.17863/CAM.68942`](https://www.repository.cam.ac.uk/items/e880c682-1876-4644-b0e5-b2c1b87cd35b)** — CC-BY-2.0-UK, 777 kB.
Authors include Tristan Youngs and Alan Soper. Contents are effectively **a complete EPSR
project**: `.mint01` files, Gudrun setup files, EPSR template molecules
(`.jmol`/`.mol`/`.ato`), a simulation box `.ato`, isotope composition
`.NWTS`/`.NWTStot.wts`, and run settings `.epsr.inp`. Raw data is not included (held at
ISIS under RB1300029).
→ **This record closes the §2.2 metadata gap for its system, and doubles as a reference
for what a well-formed EPSR project looks like.** Worth reading before designing anything.

**Dissolve examples** — [`github.com/disorderedmaterials/dissolve`](https://github.com/disorderedmaterials/dissolve),
`examples/{water,silica_bulk,cu-btc,argon,benzene,emim-oac}`. **This is where the data
already in this repo came from** (confirmed by the Dissolve-written `input.txt` headers,
the vestigial `dissolve_example_test(...)` CMake stubs, and commit `cd37a64`). The water
example is SANDALS 2001 H₂O/HDO/D₂O; silica is `SiO2_NIMROD.mint01`. The `input.txt`
files are the metadata source for composition, density and force field — i.e. a
Dissolve-input parser would be a high-value, low-risk first driver, and `argon`,
`benzene` and `emim-oac` are unused systems available immediately.

**ISIS published-data DOIs** — pattern `10.5286/ISIS.E.RB<n>`, landing pages at
`data.isis.stfc.ac.uk/doi/STUDY/<id>/`, rendering without login and showing PI,
experiment number, public release date, CC-BY-4.0 and format. **Format is "RAW/Nexus"**,
so these are raw, per Finding 2. **[CONFIRMED — landing page public; byte-level anonymous
download not exercised, UNCONFIRMED]**

**Not found:** a Soper-authored Zenodo/figshare deposit of the water structure datasets
**[UNCONFIRMED — searched, not located]**. In practice the Dissolve `examples/water`
files are the canonical public Soper-lineage water data in `.mint01` form.

### 3.2 ISIS — API confirmed working, but serves raw data

- **Catalogue:** ICAT, with DataGateway as the web front end. No SciCat migration found. **[CONFIRMED]**
- **API is live in production:** `GET https://data.isis.stfc.ac.uk/datagateway-api/ping`
  returns `DataGateway API OK`; OpenAPI spec at `/datagateway-api/openapi.json` (~1.19 MB).
  Endpoints are auto-generated per ICAT entity as `/{entity}`, `/{entity}/count`,
  `/{entity}/findone`, `/{entity}/{id}`. **[CONFIRMED by direct request]** — the specific
  paths `/investigations`, `/datafiles`, `/datasets`, `/instruments` are near-certain from
  the pattern but **[PARTIAL]** (only the first 100 kB of the spec was read).
  Source: [ral-facilities/datagateway-api](https://github.com/ral-facilities/datagateway-api) (Apache-2.0).
- **Auth:** ICAT session ID carried as `Authorization: Bearer <session_id>`, obtained from
  `POST /sessions`. No API keys, no OAuth. Missing header → 401, invalid session → 403. **[CONFIRMED]**
- **Account:** an **STFC User Office account** (username = registered email). A UKRI FedID
  works only if linked, and STFC recommends User Office credentials because some services
  reject federal logins. **[CONFIRMED]** — https://isis.analysis.stfc.ac.uk/auth/help
  Whether a non-beamtime-holder is granted an account in practice is **[UNCONFIRMED]**.
- **File downloads go through IDS** (ICAT Data Service), not DataGateway API; `python-icat`
  wraps it (`getData`, `getDataUrl`, `prepareData`) over **SOAP**, and archived data
  triggers a restore, making downloads two-phase. **[CONFIRMED]** The ISIS IDS base URL is
  **[UNCONFIRMED — not published]**.
- **Anonymous access:** ICAT ships an `anon` authenticator, but whether ISIS enables it is
  **[UNCONFIRMED]**; the PaNOSC `/search-api` returned 404/400 on the ISIS host, weak
  evidence against.
- **Embargo: three years** after the experiment ends, then released **CC-BY** with
  citation of the dataset DOI (Data Policy s.3.4.2). **Early release can be requested** —
  isisdata@stfc.ac.uk, target 72-hour response (s.3.4.6). Raw data curated ≥10 years. **[CONFIRMED]**
- **Ada** (formerly DAaaS, `data.analysis.stfc.ac.uk` → `ada.stfc.ac.uk`) is a
  browser-delivered analysis workspace with data mounted and software preinstalled —
  **not a remote data API**. Notably, a nightly GudPy build is published there, so it is
  where reduction would happen if raw→reduced ever came into scope. **[CONFIRMED]**
- **Instruments:** Disordered Materials Group (lead Daniel Bowron) operates **SANDALS**
  (hydrogenous/light elements, heavy H/D substitution, instrument DOI
  `10.5286/isis.instrument.4860`), **NIMROD** (Q = 0.01–50 Å⁻¹, interatomic to
  mesoscopic) and **GEM**; **POLARIS** is also listed for total scattering. **[CONFIRMED]**

### 3.3 ORNL (SNS/HFIR) — good discovery, no route to EPSR's format

- **ONCat API works anonymously for metadata. [CONFIRMED by direct request]**
  `GET /api/facilities` → SNS + HFIR; `GET /api/experiments?facility=SNS&instrument=NOM`
  → IPTS records with acquisition dates, run-number ranges, extensions, sizes;
  `GET /api/datafiles?...&projection=location` → paths such as
  `/SNS/NOM/IPTS-10040/0/28269/NeXus/NOM_28269_event.nxs`. Docs: https://oncat.ornl.gov/doc.
  Python client [`pyoncat`](https://pypi.org/project/pyoncat/) (GPLv3).
  → **Genuinely useful for enumerating what exists**, with no credentials.
- **But anonymous gets metadata and paths, not bytes. [CONFIRMED]** Experiment `title`,
  `users`, `members`, `tags` and `samples` are `restricted`; the Download resource states
  "Admin access required". API credentials are **issued individually** by emailing
  oncat-support@ornl.gov with your use case. OAuth scopes are documented
  (`api:read`, `api:write`, …) but the **token endpoint, grant type and client_id are
  [UNCONFIRMED — not public]**.
- **No embargo clock and no automatic release. [CONFIRMED]** Access is restricted to the
  PI and the IPTS proposal team indefinitely; "Users are responsible for the long term
  preservation of their data" (https://sns.gov/users/data-management). Beam time is free
  in exchange for *publishing*, not for releasing data. Public ORNL neutron data exists
  only where a PI opted in to mint a DOI (prefixes `10.13139/ORNLNCCS/*`,
  `10.14461/oncat.data/*`). **No published DOE/ORNL neutron embargo period was found
  [UNCONFIRMED — likely none exists]**.
- **Instruments:** **NOMAD** (BL-1B, purpose-built for PDF, typical Q_max 30–35 Å⁻¹) and
  **POWGEN**. **[CONFIRMED]**
- **Reduction chain and output:** event NeXus → Mantid + `mantid_total_scattering` →
  ADDIE → pystog. Outputs: `SofQ/` and `SofQ_merged/` as **NeXus**, plus plain-text
  `StoG/*.sq`, `*.gr`, `*_rmc.fq`, `*_rmc.gr` (RMCProfile conventions). **[CONFIRMED]**

### 3.4 Why ORNL is harder than ISIS — and it is a physics problem, not a parsing one

ISIS total-scattering data reduces through **Gudrun** (Soper, RAL-TR-2011-013), the same
lineage as EPSR, and **`.mint01` is literally Gudrun's output** — the merged differential
cross-section after Top-Hat post-processing and real-space deconvolution below `r_min`.
Gudrun ships predefined settings for **SANDALS, NIMROD, GEM, POLARIS (ISIS) and D4c
(ILL)** — **no SNS or HFIR instrument**. **[CONFIRMED]**

Converting ORNL output would be easy as text and risky as physics:

- ORNL `SofQ` is **normalised S(Q)**, documented as approaching 1 at high Q. **[CONFIRMED]**
- Gudrun `.mint01` is **not normalised** — corroborated inside this repo by Dissolve's own
  import directive for these very files, `ReferenceNormalisedTo None`
  ([`experimental_data/water/input.txt:158-166`](experimental_data/water/input.txt#L158-L166)),
  and by the negative low-Q values in every `.mint01` here. **[CONFIRMED]**
- RMCProfile's own guide warns "F(Q)" is used loosely in the community — some mean
  Q[S(Q)−1] — and advises checking the low-Q region to identify the convention. The same
  ambiguity bites an ORNL→EPSR conversion. **[CONFIRMED]**
- **No official ORNL→Gudrun/EPSR converter and no documented ORNL→`.mint01` path was
  found [UNCONFIRMED — believed not to exist]**. ORNL's analysis-software page lists EPSR
  and Dissolve as downstream tools but specifies no input formats and never mentions Gudrun.

→ **Treat ORNL as secondary.** If implemented, the conversion needs a written statement of
the assumed convention and a physics review — ideally validated by reducing one ISIS
dataset both ways and comparing.

### 3.5 Other facilities

- **ILL** — the strongest second facility, because **Gudrun natively supports D4c**.
  Three-year non-disclosure, **extended to five years where no request for the data has
  been made**, then public via the ILL Data Portal (data.ill.eu) with DOIs
  `10.5291/ILL-DATA.*`. A programmatic API is **[UNCONFIRMED]**. **[CONFIRMED policy]**
- **PSI (SINQ)** — data policy **[UNCONFIRMED]**; PSI originated **SciCat**, so a
  SciCat API likely exists **[UNCONFIRMED]**.
- **J-PARC MLF** — **NOVA (BL21)** is the total-scattering instrument; policy and remote
  access **[UNCONFIRMED]**.

---

## 4. Proposed architecture

Four stages, so each is independently testable and the credentialed parts are isolated.

```
  discovery ──► retrieval ──► conversion ──► validation ──► (metadata resolver) ──► run folder
   search()      fetch()      to Q,S,E        contract §2.1    human-reviewable
```

**`DataSource` driver interface** — one implementation per repository:

```
search(query)      -> [Record]      # Record: id, title, instrument, doi, licence, files, metadata
fetch(record, dir) -> [LocalFile]   # downloads; returns paths + checksums
describe(record)   -> Metadata      # whatever composition/isotope/FF info the source carries
capabilities()     -> {needs_auth, provides_reduced, provides_metadata, licence_known}
```

`capabilities()` is what keeps Finding 2 from being rediscovered the hard way: an ISIS
catalogue driver declares `provides_reduced = false`, and the agent can say so up front
instead of fetching gigabytes of raw NeXus that EPSR cannot read.

**Suggested drivers, in build order:** `DissolveExamples` (git/HTTP, no auth, ships
metadata) → `DataciteArchive` (generic DOI-resolving driver covering Bath/Cambridge/Zenodo,
no auth) → `LocalDirectory` (private/unpublished data a user already has — satisfies the
"private sources" requirement with no credential handling at all) → `IsisIcat` (optional
credentials, discovery + raw) → `OrnlOncat` (anonymous discovery; gated retrieval).

**Auth provider** — optional per driver, resolved from environment or a system keyring,
**never from code, never committed, never logged**. A driver that needs credentials and
has none must degrade to "unavailable" with a clear message, not fail obscurely.

**Converter** — normalises any supported input to `Q, S(Q), E(Q)`. Must record the assumed
normalisation convention in the output header as a `#` comment (free, and prevents
silent convention mistakes later).

**Validator** — enforces §2.1 mechanically: header-token rule, three parseable columns,
Q monotonic ascending, ≤5000 points, blank-free filename, trailing sentinel row present,
and a **warning** when data extends beyond `nq × qstep` (silently discarded) or stops
short of it. The four `.mint01` files already in the repo are ready-made fixtures.

**Metadata resolver** — explicitly separate and human-in-the-loop: ingest metadata where
the source provides it (Dissolve `input.txt`, the Cambridge record's `.mol`/`.NWTS`),
otherwise present a candidate for confirmation. It should emit the `.mol`, `mixato`
answers and `.NWTS.dat`, then drive the existing chain in
[`validation/PYTHON_WRAPPER_WORKFLOW.md`](validation/PYTHON_WRAPPER_WORKFLOW.md) §1
rather than writing `.ato` files directly.

---

## 5. Provenance — a cheap, high-value addition

**Today nothing records where any dataset came from.** There is no manifest, no checksum
file, no DOI and no retrieval date for any of the four `.mint01` files in this repo;
`rust/docs/case_summaries/references.json` lists papers but no dataset URLs. Origin is
recoverable only by reading Gudrun's own header lines (run numbers, experimenter,
reduction date) and inferring the rest from git commit messages.

The agent should write a manifest beside each acquired dataset: **source driver, record
ID/DOI, instrument and run numbers, retrieval timestamp, checksum, licence, and the
conversion applied.** That makes new validation cases as auditable as the existing ones
and is nearly free to implement. Retro-fitting manifests for the four existing files is a
sensible first task — the information is still recoverable from the file headers.

---

## 6. Constraints and risks

- **Credentials**: never committed, never logged, env/keyring only; each driver works or
  cleanly declares itself unavailable.
- **Terms of use**: use the documented APIs, not portal scraping. ISIS access is
  "restricted to those who register with STFC/ISIS as users" (Data Policy s.3.6.1); ORNL
  asks you to describe your integration before issuing API credentials. Honour rate limits.
- **Embargo**: ISIS data is restricted to the experimental team for three years; ORNL data
  is restricted indefinitely unless a DOI was minted. An agent must not attempt to reach
  around this — the supported route for early ISIS access is an email request.
- **Licensing for redistribution**: post-embargo ISIS data is CC-BY; the Bath record is
  CC-BY-4.0; the Cambridge record is CC-BY-2.0-UK — all require attribution if committed
  into this repo. Note separately that the **EPSR distribution's own terms state
  "Demonstration and raw data files are for testing only"** and commercial use needs
  permission ([github.com/disorderedmaterials/EPSR](https://github.com/disorderedmaterials/EPSR),
  no standard licence) — relevant to `EPSR-test/` and `test_problems/`, not to data the
  agent fetches.
- **Silent data loss**: Q beyond `nq × qstep` (≈30 Å⁻¹ in the validation cases) is read and
  discarded without warning; the `.mint01` files here extend to ~50 Å⁻¹. The validator
  should surface this.
- **Operational quirks the agent must respect** (all documented in
  [`validation/PYTHON_WRAPPER_WORKFLOW.md`](validation/PYTHON_WRAPPER_WORKFLOW.md) §6):
  the setup tools read/write `system_commands.txt`/`dirlist.txt` **in the current working
  directory**, so a driver must `chdir` into the case folder and parallel setups need
  separate processes; `fmole` with the template `.mol` still present is unsupported in the
  Rust port; start with **no** `.pcof` rather than an empty one; and never regenerate the
  `.EPSR.inp` between runs, because EPSR writes `ireset 0`/`iinit 0` back into it.

---

## 7. Recommended build order

1. **Converter + validator** against the four existing `.mint01` files as fixtures.
   No network, no credentials, immediately testable. Reuse the reader in
   `EPSR-test-ex95/scripts/make_data.py`.
2. **`DissolveExamples` driver** — proves the interface end-to-end with zero auth, and
   immediately yields unused systems (`argon`, `benzene`, `emim-oac`, `cu-btc`) plus the
   `input.txt` metadata the resolver needs.
3. **Manifest/provenance layer**, retro-fitted to the existing four files.
4. **`DataciteArchive` driver** — start with Cambridge `10.17863/CAM.68942`, since it
   contains a full EPSR project and therefore also validates the metadata resolver, then
   Bath `10.15125/BATH-00359` for a multi-sample H/D series.
5. **`LocalDirectory` driver** — covers private/unpublished data with no credential path.
6. **`IsisIcat` driver** — optional credentials, `provides_reduced = false`; valuable for
   discovery, provenance and DOI resolution, and the foundation for raw→reduced later.
7. **`OrnlOncat` driver** — anonymous discovery first (it already works); retrieval and
   the S(Q) convention conversion only after a documented physics review.

Steps 1–5 need no facility account at all, which matches the "works without credentials"
requirement and delivers most of the practical value.

### If raw→reduced is ever brought into scope

That is the step that unlocks ISIS's actual holdings. It would mean driving **Gudrun**
([github.com/disorderedmaterials/Gudrun](https://github.com/disorderedmaterials/Gudrun),
Fortran, CMake/HDF5; no stated licence) or **GudPy**
([github.com/disorderedmaterials/GudPy](https://github.com/disorderedmaterials/GudPy),
whose "Export Archive" emits exactly a zip of `.mint01` files) — and supplying per-experiment
expert choices: detector purging, attenuation and multiple-scattering corrections, deadtime,
monitor normalisation, background and container subtraction, **vanadium** absolute-scale
normalisation, **Placzek**/inelasticity modelling, self-scattering subtraction, then Top-Hat
smoothing and real-space deconvolution (RAL-TR-2011-013 ch. 2–3). It is a substantial
project in its own right and should not be bolted onto the acquisition agent.
