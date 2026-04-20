# Stochastic Microtubule Dynamics with EB1 Tip Tracking

Previous to MT-archProt

A Monte Carlo Gillespie-based stochastic simulation of microtubule dynamics with MT-bound proteins, including GTP hydrolysis and monomeric protein (EB1-like) protein binding and unbinding.
Based on 

---

## File Structure

| File | Purpose |
|---|---|
| `params.py` | All simulation parameters and constants |
| `initialization.py` | Builds all data structures before the main loop |
| `helpers.py` | Geometry utilities and pocket classification |
| `simulation_helpers.py` | State query and update functions |
| `simulation.py` | Main Gillespie loop, event generation and execution |

---

## How to Run

```bash
python simulation.py
```

Output is a `simulation.gif` saved in the working directory.

---

## Parameters

All parameters are defined in `params.py`, including:

| Parameter | Description | Default |
|---|---|---|
| `n_pf` | Number of protofilaments | 13 |
| `n_iterations` | Total Gillespie steps | 450000 |
| `seed_length` | Starting seed length in dimers | 25 |
| `conc_tubGTP` | Free GTP-tubulin concentration (uM) | 12.0 |
| `conc_prot` | EB1 concentration (uM) | 0.200 |
| `kon_tub` | Tubulin on-rate per PF (uM-1 s-1) | 0.65 |
| `koff_tubGTP` | Tubulin off-rate, GTP below (s-1) | 0.2 |
| `koff_tubGDP` | Tubulin off-rate, GDP below (s-1) | 200.0 |
| `k_hydrolysis` | GTP hydrolysis rate per PF (s-1) | 0.0423 |
| `k_lateralbond` | Lateral bond formation rate (s-1) | 100.0 |
| `taper_threshold` | Max taper before catastrophe trigger (dimers) | 75 |

See `params.py` for the full list including all lateral bond breakage rates and EB1 on/off rates.

---

## Data Structures

### `MT_lattice` — shape `(13, 300, 2)`
Tracks the state of every tubulin dimer.
- `axis 0` — protofilament index (0-based)
- `axis 1` — height position (0 = minus-end base)
- `axis 2, layer 0` — hydrolysis state: `0` = GTP, `1` = GDP
- `axis 2, layer 1` — right lateral bond count: `0`, `1`, or `2`

### `prot_sites` — shape `(14, 300, 3)`
Tracks the state of every EB1 binding pocket (grooves between protofilaments).
- `axis 0` — groove index (groove `g` sits between PF `g-1` and PF `g`)
- `axis 1` — height position
- `axis 2, layer 0` — site type: `SITE_LATTICE`, `SITE_EDGELAT2`, `SITE_EDGELONG2`, `SITE_EDGE3`, `SITE_SINGLE`, `SITE_SEAM`
- `axis 2, layer 1` — nucleotide state: `NUC_GTP`, `NUC_MIXED`, `NUC_GDP`
- `axis 2, layer 2` — EB1 bound: `0` = unoccupied, `1` = occupied

### `pf_len` — shape `(13,)`
Length of each protofilament in dimers.

### `highest_lat` — shape `(13,)`
For each PF, the height of the highest subunit that is right-laterally bonded.

### `highest_full_GDP` — integer
The highest height at which all 13 protofilaments are GDP-hydrolyzed. Used to limit the inner hydrolysis loop.

---

## Events

Each Gillespie step generates a pool of candidate events and picks the one with the shortest sampled waiting time. Events are:

| Event | Description |
|---|---|
| `tub_add` | Tubulin addition at PF tip |
| `tub_removal` | Tubulin dissociation from PF tip |
| `lat_bond_form` | Lateral bond formation |
| `lat_bond_break` | Lateral bond breakage |
| `prot_bind_edge` | EB1 binding at an edge site |
| `prot_bind_lattice` | EB1 binding at a lattice site |
| `prot_remove_gtp_edge` | EB1 unbinding from GTP edge site |
| `prot_remove_gdp_edge` | EB1 unbinding from GDP edge site |
| `prot_remove_gtp_lattice` | EB1 unbinding from GTP lattice site |
| `prot_remove_gdp_lattice` | EB1 unbinding from GDP lattice site |

### GTP Hydrolysis
Hydrolysis is not a Gillespie event. Instead, after each winning event fires, every non-terminal GTP dimer above the seed is checked. A waiting time is sampled from `Exp(k_hydrolysis)` and if it is less than the winning event's `dt`, the dimer hydrolyzes. This is computationally more efficient than adding one event per GTP dimer to the candidate pool.

---

## Output

Every `snapshot_freq` steps the simulation records:
- `out_time` — elapsed simulation time
- `out_pf_lengths` — full PF length array
- `out_n_bound` — total EB1 proteins bound

A GIF is saved showing two panels at each snapshot:
- **Left** — bar chart of protofilament lengths
- **Right** — lattice map with GTP (medium blue), GDP (dark navy), EB1 proteins (green dots), lateral bonds (orange lines), and seam bonds (yellow/red lines)
