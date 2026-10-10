# Opportunistic network-level ISAC: reproduction and model audit

A runnable terrestrial-model study of **“Opportunistic Network-Level ISAC with Cooperative Sensing: A Meta-Distribution Analysis”**, Yasser Nabil, Hesham ElSawy, and Hossam S. Hassanein, [arXiv:2604.06069v1](https://arxiv.org/abs/2604.06069).

**Status: preliminary reproduction study with validated numerical components. No claim of exact reproduction.** The full primary paper was read before model implementation. The necessary communication derivations in [reference 3, arXiv:2502.13348v2](https://arxiv.org/abs/2502.13348) were also read. Original Figure 3 data are unavailable. Important printed-equation/conditioning discrepancies and unspecified settings remain; see the report and audit.

Start reading:

- [Beginner tutorial, notation, SI parameters, and equation derivations](docs/tutorial.md)
- [Model audit and unresolved questions](docs/model_audit.md)
- [Independent monostatic joint-interference derivation](docs/monostatic_derivation.md)
- [Results and validation report](docs/report.md)
- [Source versions and PDF fingerprints](docs/sources.json)

## Windows PowerShell

Python **3.12** and the dependencies below were verified in the cloud environment. No Unix-only shell commands are needed for the Python workflow.

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m src.isac_reproduction.check_environment

# Small, runnable monostatic baseline first
.\.venv\Scripts\python.exe -m src.isac_reproduction.cli baseline --samples 1000 --qmc-power 8 --output results/baseline

# Coverage: original paper's 200000 network realizations / 25 km²
.\.venv\Scripts\python.exe -m src.isac_reproduction.cli coverage --samples 200000 --qmc-power 8

# Meta-distribution: explicit distance conditioning and a marked-geometry comparison
.\.venv\Scripts\python.exe -m src.isac_reproduction.cli meta --outer 128 --slots 128 --qmc-power 9

# Density and energy-split sweeps; sample count is preliminary
.\.venv\Scripts\python.exe -m src.isac_reproduction.cli rates --samples 6000

# Area, sample-count, spatial-integration, numerical and model checks
.\.venv\Scripts\python.exe -m src.isac_reproduction.cli convergence --samples 10000
.\.venv\Scripts\python.exe -m src.isac_reproduction.cli validate

# Focused inversion/transform and positions-only conditioning follow-ups
.\.venv\Scripts\python.exe -m src.isac_reproduction.inspection
.\.venv\Scripts\python.exe -m src.isac_reproduction.positions_meta
```

On Linux replace `.\.venv\Scripts\python.exe` with `.venv/bin/python` after `python -m venv .venv`. Running with the already installed scientific Python also works: `python -m src.isac_reproduction.cli baseline --samples 1000`.

No special credentials or external services are needed for simulation. The plotting module configures writable temporary caches and a headless renderer. It saves PNG and SVG artifacts; no display server is needed. The independent dependency check can be run with `MPLCONFIGDIR` and `XDG_CACHE_HOME` set to writable directories if your home is read-only.

All commands run from the repository root. `--help` lists options. Use `--output` to preserve a previous run rather than overwrite its files. Outputs contain configuration, seed, actual counts, diagnostic bounds, and source-code hashes. Binary NPZ files contain numerical samples, not just figures. All RNGs are explicitly seeded; changed sample counts may change the random-stream sequence. Statistical agreement rather than bitwise identity across SciPy/platform versions is the relevant check.

`requirements.txt` pins the three requested scientific libraries. For the complete verified dependency set, install `requirements-lock.txt` instead. It records transitive dependency versions from the fresh Python 3.12 environment. The two focused follow-up commands expect the default-model `meta` artifacts in `results/`; they retain preliminary nested-MC labels.

## Configuration and model choices

The centralized [Config](src/isac_reproduction/config.py) contains all §IV values in SI units. [paper_defaults.json](configs/paper_defaults.json) lists configurable input fields. Partial JSON overrides are supported:

```powershell
.\.venv\Scripts\python.exe -m src.isac_reproduction.cli baseline --config configs/per_bs_budget.json --output results/per_bs
```

The default interprets the 1 W budget **per beam**; the alternative splits it among eight beams. This scope is not specified unambiguously by the paper. Uniformly rotated equally spaced non-overlapping beams are another disclosed implementation choice. Serving target/user links are forced LoS, not selected by nearest-LoS association. Sensing uses the paper's stated top-hat simplification; communication uses its cosine pattern. No NTN extension, path-loss regularization, hidden noise figure, or Rmax distance filtering is introduced.

Simulation draws independent exponential original echo RCS powers. It never substitutes the analytical Gamma approximation. It retains common point locations and beam rotations for direct and reflected paths, and retains target marks when echoes move from interference to useful power as cluster size increases.

The printed analytical model separately implements Eqs. (5),(6),(11)–(20). Its nearest-point omission and incomplete distance conditioning remain **visible**; the independent monostatic audit is a separately labeled derivation. Communication implements [3]'s cosine-pattern Alzer derivation in both frequency-separated and shared-band interpretations because the primary paper changes the frequency reuse architecture. Its chosen rate curves use the primary paper's shared-band statement.

## Outputs and limits

| Artifact | What it contains |
|---|---|
| `results/figure3a.{png,svg,npz,json}` | Original exponential-echo simulation, 95% coverage intervals, printed analysis, independent monostatic audit |
| `results/figure3b.{png,svg,json}`, `meta_nc*.npz` | Distance-conditioned reliability, raw Gil–Pelaez inversion checks, nested samples |
| `results/conditioning_comparison.{png,svg}` | Distances versus full **marked** geometry; different conditioning |
| `results/figure3c.{png,svg}`, `rates_density.npz` | Density sweep, all five sensing cluster sizes and both communication interpretations |
| `results/figure3d.{png,svg}`, `rates_gamma.npz` | Energy split, sensing/communication sums and communication-only comparison |
| `results/convergence.json`, `area_*.npz` | 6.25/25/100 km² sensitivity, quarter/full sample precision, 64/256/1024 spatial samples |
| `results/validation.json` | Echo moments, Gamma/Alzer errors, transforms, series, inversion, limiting cases, communication checks, joint monostatic audit |

Figure 3(b)'s nested estimate is deliberately **preliminary**: 128 inner samples can give a worst-case 95% half-width of about 0.087 for a single conditional probability. An empirical meta-CCDF of these estimates inherits that uncertainty, especially near reliability one. Use `--outer 512 --slots 1024` for a more expensive follow-up; do not confuse rounded inner success fractions with exact reliability-one atoms. The plotted analytical CCDF is computed from conditional probabilities; raw inversion output is retained without clipping.

Coverage and rate confidence intervals quantify sampling uncertainty, not modeling ambiguity. Direct/reflected LT validation against their own surrogate processes does not validate factorization for a shared physical PPP. Original-paper numerical agreement remains unresolved. Read [the report](docs/report.md) before using the figures in a publication.
