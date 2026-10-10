# Opportunistic network-level ISAC reproduction

Status: **project scaffold; paper access blocked; no model implemented or validated**.

Target: “Opportunistic Network-Level ISAC with Cooperative Sensing: A Meta-Distribution Analysis”, Yasser Nabil, Hesham ElSawy, and Hossam S. Hassanein, https://arxiv.org/pdf/2604.06069.

The full paper must be read before implementing its terrestrial model. No paper parameters, equations, or numerical results are inferred in this scaffold. Reference [3] must be identified from the supplied paper and consulted where necessary.

## Windows PowerShell

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m src.isac_reproduction.check_environment
```

On Linux use `python -m venv .venv`, `.venv/bin/python -m pip install -r requirements.txt`, and `.venv/bin/python -m src.isac_reproduction.check_environment`.

This command checks only the scientific Python dependencies; it is not a simulation or reproduction check.

## Planned validation

1. Read full paper; extract equation references, units, parameters, and unresolved assumptions.
2. Implement original random monostatic echo powers, then cooperative sensing with shared geometry.
3. Implement analytical approximations independently; distinguish distance conditioning from full geometry conditioning.
4. Validate interference transforms, moments, coverage, rates, and inversion, including finite-area and sample-count sensitivity.
5. Save seeded numerical results and confidence intervals before plotting Figure 3(a–d).

Use SI units internally. Never replace missing paper settings with silent defaults. Results will be labeled preliminary, validated, or reproduced according to the evidence.
