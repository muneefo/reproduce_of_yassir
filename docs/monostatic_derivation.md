# An independent monostatic audit

This derivation is an additional diagnostic, **not a replacement attributed to the paper**. It is exact for an infinite PPP with the implemented top-hat sensing beams, equally spaced rotations, independent link-blockage marks, and the paper's exponential RCS powers. Numerical integration and spatial averaging still require convergence checks.

## Why derive it?

Nc=1 needs no Gamma moment approximation and no Alzer approximation: its desired echo is exponential with mean μ1. Given R1, the target-centered association disk is empty of other BSs. An exact conditional survival probability can therefore be obtained from a **joint** interference PGFL. It tests whether a discrepancy in the printed analysis could be blamed on Gamma matching. For monostatic sensing it cannot.

## Geometry first

Put the target at (0,0) and its sensing BS at (R1,0). The receiver points left toward the target. For a direct interferer in its receive sector, write position (R1−d cosφ,−d sinφ), with |φ|≤θb/2. Its distance from the target is R=√(R1²+d²−2R1 d cosφ).

The point is outside the empty target disk precisely when d≥2R1 cosφ. Thus direct interference has a **sector-dependent exclusion**, rather than the independently sampled circular guard in Eq. (11). This follows directly by setting R²≥R1² and solving d²−2R1d cosφ≥0.

## Keep the shared beam marks

The illumination marginal of the equally spaced random beam set is pB=Nbθb/(2π). Let δ be the smallest periodic angular separation, modulo 2π/Nb, between the directions from the interferer to the target and to the sensing BS. The probability both directions are illuminated is

    pBoth = max(θb − |δ|, 0) / (2π/Nb).

The max operation here measures the length of intersection of two angular intervals; it is not clipping a numerical probability result. This formula uses θb≤π/Nb, as satisfied by the stated parameters.

For Laplace argument s, let a=sAσav,bi/R1². After averaging independent target RCS/blockage, define

    fT(R) = exp(−ηb R) · a/(R²+a).

After averaging independent direct fading/blockage, define

    fD(d) = pLOS(d)[1−(1+sPcGm²CL d^−αL/mL)^−mL]
            + (1−pLOS(d))[1−(1+sPcGm²CN d^−αN/mN)^−mN].

The per-point joint LT loss inside the receiver sector is

    pB fT + pB fD − pBoth fT fD.

This comes from expanding E[(1−Btarget fT)(1−Btx fD)]. It explicitly retains the covariance between the two illumination indicators. Outside the receiver sector only pB fT remains.

## Apply PGFL and the exponential desired tail

The reflected-only exponent over the entire exterior of the target disk is

    ET = λBS Nbθb ∫_(R1)^∞ exp(−ηb r) a r/(r²+a) dr.

It starts at R1 for monostatic sensing and **naturally includes** the first outside PPP point. The extra direct/joint contribution is

    ED = λBS ∫_(−θb/2)^(θb/2) ∫_(2R1cosφ)^∞
           [pB fD(d) − pBoth fT(R(d,φ)) fD(d)] d dd dφ.

The desired exponential power gives s=τ/μ1. Therefore

    P(SINR>τ | R1) = exp[−ET − ED − s(KBT W + PcζSI)].

`monostatic_audit.conditional_coverage` evaluates this expression using the validated E1 integral for ET and Gauss–Legendre quadrature for ED, mapping the radial interval to [0,1]. The infinite-area PPP is used; the Monte Carlo remains finite-area and is checked for sensitivity. The numerical source and approximation labels are explicit in plots and JSON reports.

This audit isolates a useful fact: monostatic success can agree with the geometry simulation while Eq. (14), implemented literally using (11) and (13), does not. It does not identify which undocumented implementation the authors used to generate their plotted analytical curves. That requires their numerical data/code or clarification.
