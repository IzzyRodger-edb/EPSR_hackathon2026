# Experimental partial pair distribution functions of water

Downloaded on 2026-10-17 from the ISIS Disordered Materials Neutron Database:
http://wwwisis2.isis.rl.ac.uk/Disordered/Database/DBMain.htm → Water Data
(`New Water/DB_NewH2O.htm`).

| file | content |
|---|---|
| `h2o298aoogor.txt` | g_OO(r), H2O at 298 K, case (a) |
| `h2o298aohgor.txt` | g_OH(r), H2O at 298 K, case (a) |
| `h2o298ahhgor.txt` | g_HH(r), H2O at 298 K, case (a) |

Reference: A. K. Soper, "The radial distribution functions of water and ice from
220 to 673 K and at pressures up to 400 MPa", *Chem. Phys.* **258** (2000) 121–137.

Columns: bin number, r (Å), g(r), standard deviation. Each file ends with an
all-zero terminator row, which `epsrvis.read_isis_gr` drops (also `epsr.vis.read_isis_gr` in the wrapper). The partials were
extracted from H₂O/D₂O/HDO neutron diffraction data by an EPSR-type
refinement, so they are model-assisted, not directly measured. The partials include the
intramolecular O–H and H–H peaks.

**Why case (a).** The paper gives two cases at 298 K. Case (a) has the
ambient-water O–O structure: first peak at 2.73 Å (height 2.75), minimum at
3.36 Å, second peak at 4.5 Å. Case (b) has its second O–O shell pulled in to
4.26 Å, the signature of compressed water. Case (a) is therefore the ambient
pressure state, which matches the EPSR simulation (300 K, 0.1 atoms/Å³). The
database page does not label the cases, so this assignment comes from the data.

**Silica.** The database has no experimental *partials* for SiO₂, only the total
correlation function T(r) of vitreous silica (`sio2tor.txt`). So the silica
assets have no experimental partials.
