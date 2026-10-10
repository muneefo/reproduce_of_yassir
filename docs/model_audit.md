# Model audit: definitions, approximations, and choices

Read this before interpreting a curve as reproduction. The implementation keeps the printed analytical equations separate from the shared-geometry Monte Carlo. Differences are evidence to investigate, not permission to retune parameters.

## Conditioning

The paper explicitly defines PS(R)=P(SINR>τ | R1,…,RNc). This integrates out outside points, cluster azimuths, beam orientations, blockages, fading, and RCS marks. It is **not** a success probability conditional on the whole PPP configuration, and it is also not the usual fading-only reliability conditional on persistent beam/blockage marks. `conditional_simulation(mode="distances")` fixes ordered radii and resamples the other quantities. `mode="full"` fixes an entire point pattern, rotations, and blockages, and resamples fading/RCS; it is a distinct, explicitly marked comparison.

`mode="positions"` conditions on all BS positions alone and resamples all beam/blockage/fading/RCS marks. This is the complete **unmarked PPP geometry** conditioning, distinct from both ordered-distance conditioning and complete marked-geometry conditioning. `positions_meta.py` evaluates it and saves `positions_conditioning.png` and per-cluster numerical samples. The paper's written definition remains the ordered-distance one.

The law of total expectation gives E[PS(R)]=P(SINR>τ). For the same sample, averaging group success fractions equals pooling group successes exactly. That identity is an implementation consistency check, **not independent statistical validation**. Compare against the independently seeded unconditioned `figure3a` estimate and use the outer-group confidence interval. Different meta-distributions can have the same first moment because conditioning changes dispersion.

## Nearest interferers and two exclusion centers

1. Target-centered nearest association gives an empty disk B(target,R1); it does **not** give a receiver-centered empty disk of radius R1. The sensing BS is displaced from the target. Conditional cluster locations contribute direct interference too. Equation (11) replaces this geometry with an independent receiver-centered guard RI, whose PDF uses unthinned λBS. If RI means nearest **aligned** interferer, its PDF would use λI, not λBS. The printed PDF is instead the unthinned nearest-potential-BS law.
2. If RI is the nearest potential interferer, that BS is a point on the random boundary. PGFL strictly outside RI does not include its fading/LoS/alignment contribution. Equation (11) has no boundary-point factor. If RI is instead an independently imposed random guard, (11) is a legitimate surrogate LT; it is not the exact receiver-centered nearest-interferer model. We validate its integral and surrogate Monte Carlo separately.
3. Equation (8) includes Ri≥R_(Nc+1), so the (Nc+1)-st BS should contribute when illuminated. Equation (13) uses only a PPP integral starting at that radius. A conditioned nearest point is not supplied by a measure-zero integration boundary; its Bernoulli/exponential mark requires an explicit factor. The literal implementation therefore omits its contribution. For a true PPP conditioned on the full ordered cluster R, an equivalent clean reflected-interference LT is PGFL outside **RNc**, with no separate boundary point. It includes the first subsequent PPP point naturally.
4. Equation (12) is correct for R_(Nc+1)|R1. Equation (14) conditions on **all R**, but its reflected transform still averages (12) rather than the appropriate R_(Nc+1)|RNc. Thus its reflected interference does not depend on RNc except through cluster size. This is not the same conditional model. The code retains the printed transform; `validation.json` contrasts it with a fixed-R1 exact reflected PPP baseline.

## Correlations and independence

The original point pattern is shared by useful echoes, direct interference, and reflected interference. Nonserving beam rotations are shared for directions toward the target and toward the sensing BS. Target-link LoS and RCS marks are reused when a cluster-size change reclassifies an echo. Cluster beams/RCS that raise desired power can also align with direct interference. Increasing Nc does not remove those BSs from direct interference.

Conditional on positions, direct-link and target-link blockage marks are independent draws: no spatial blocker process is specified. This is a transparent modeling choice; it does not establish physical blockage independence. Echo RCS marks are independent across BSs/slots as specified. Full-model useful power and network interference need not be independent after averaging beam geometry. LoS/NLoS thinning of one PPP produces independent PPPs unconditionally under independent marks, but mixing over a shared random guard couples their marginal LTs. Equation (11) is implemented as its joint guard average, not as a product of separately averaged LoS and NLoS LTs.

Equation (14) multiplies the printed reflected and direct transforms, assuming their independence and separation from the useful-power approximation. Simulation does not make this replacement. A small Pearson correlation would not prove independence; rare strong events and nonlinear Laplace expectations matter. We record observed correlations but do not treat them as validation of factorization.

## Antennas and energy

The full cosine support is 2θb, while the sensing top-hat support is θb. The primary paper states the top-hat sensing approximation in the paragraph following (4). We use that sensing pattern in **both** the geometry simulation and printed sensing analysis, and the cosine in communication. Whether the authors' Fig. 3 simulation used the full cosine for sensing is not specified. Therefore our simulation is geometry-preserving for the stated sensing powers/top-hat model, not a claim to recover every unspecified antenna implementation detail of the authors' simulator.

For equally spaced Nb=8 beams, boresights are 45° apart; 18° sensing top-hats and 36° cosine supports do not overlap. This is one concrete rotationally invariant layout with the correct illumination marginal in (4). The paper only specifies non-overlap, not the joint beam-spacing distribution. Directional correlations can depend on that layout. The sensing receiver is steered at the target. Nominal Gm and θb are held fixed during density/γ sweeps; no gain-width tradeoff is introduced.

The 1 W slot-averaged budget is not explicitly labeled per beam versus per BS. Default `power_scope="per_beam"` applies the §II powers directly in each beam equation. `"per_bs"` divides E among Nb beams. This changes noise-relative coverage even where interference ratios mostly cancel. Residual SI is Pc ζSI as printed, without extra multiplication by Nb; cross-beam SI is not specified. No hidden noise figure, receiver loss, pathloss cutoff, extra Nakagami echo factor, or Rmax rejection is introduced.

The selected exact physical constants differ slightly from possible rounded c=3×10⁸ conventions. This difference is minor but disclosed. Rmax sets Tt, not a hard cutoff on useful/interfering distances in Eqs. (8),(15),(16). The nearest distance exceeds Rmax with probability exp(−9)≈0.0001234. Practical timing rejection of such links would be another model, and is not silently imposed.

## Gamma fitting and Alzer

The useful-power sum is not generally Gamma. A Bernoulli-exponential bistatic term even has an atom at zero. Our simulated power samples use the original random variables. Moment fitting and Alzer survival substitution are isolated in `analysis.py`.

For the stated means, μbi,n/μmono=0.7(R1/Rn)^αL≤0.7. Each bistatic term contributes at least 2pn μn[μmono−(1−pn)μn]≥0 to mean²−variance, with additional nonnegative cross terms, so **keff≥1** for the stated parameters. This explains why the paper's chosen Alzer regime applies here. Arbitrary altered RCS settings need not preserve that regime. Noninteger keff must not be rounded, truncated at ceil(keff), or fed to an integer-only factorial sum. The generalized series is infinite and converges slowly near zero total interference. We use a recurrence and a conservative tail bound, report the achieved bound, and test tighter tolerances. The default numerical tolerance is 2×10⁻⁴ per conditional probability, not an estimate of model-approximation error.

Alzer with the specified coefficient is a bound-related survival approximation, not an exact Gamma CDF outside k=1. For k>1 it tends to overestimate survival. `validation.json` separates original-echo-versus-Gamma error from Gamma-versus-Alzer error. Gil–Pelaez finite frequency truncation can ring, especially near empirical jumps. We do not clip those values or call them exact inversion. Complex moments are evaluated only for finite positive conditional probabilities; actual zeros require explicit atom handling, not a fake epsilon.

## Communication transfer from [3]

Retrieved https://arxiv.org/pdf/2502.13348 (v2). Read §II-B, §IV, and Appendices F–G. [3] assigns distinct frequency blocks to beams. Thus an interfering BS has one cochannel beam; its Lemma 7 has coefficient λBS times the integral over cosine support. The primary paper instead states all Nb beams share W. For a random equally spaced non-overlapping beam set, the same-band orientation integral has coefficient **Nb λBS**. These alternatives are `ref3_single` and `same_band`. The latter follows the primary paper's physical SDMA statement; the former is the literal interference multiplicity from the cited derivation. Both are evaluated and labeled. No misalignment distribution is imported: Eq. (21) explicitly gives desired Gm. Clutter, different-frequency receive fusion, bistatic-angle RCS scaling, and selection combining from [3] are not part of the primary paper's new power-combining model.

## Finite-area and inference choices

The target is at the center of a finite square; points are Poisson with mean λBS A and uniform in that square. No wraparound introduces artificial copies. The nearest BS can be off-center but normally far from the square edge. All simulated components share this geometry. Conditioned distance simulation inserts the ordered cluster at random azimuths and samples the PPP outside its target-centered disk, within the square. A cluster extending beyond the inscribed circle raises an error instead of censoring it.

Remote direct NLoS tails decay as ∫r^(1−αN)dr, and LoS/reflected tails also have exponential blockage attenuation. We compare 6.25,25,100 km² to test truncation empirically rather than assume it negligible. Area differences within confidence intervals do not prove zero truncation error. The saved sample-count prefixes show Monte Carlo precision; spatial QMC, radial integration, guard integration, series truncation, and inversion truncation are tested separately.

Coverage Wilson intervals use independent unconditioned network realizations. Rate intervals are Student mean intervals over log1p(SINR). Conditional-distance mean intervals use independent **outer groups**. Nested meta-CCDFs have inner binomial-estimation error, particularly near t=1, and do not have the precision of exact PS(R) samples. Requiring 200000 realizations does not remove analytical conditioning mismatches.

## What remains unresolved with the authors

- Confirm per-beam versus per-BS 1 W budget, receive/cross-beam SI treatment, beam placement law, and sensing simulation antenna pattern.
- Clarify nearest-interferer boundary factors and how conditioning on RNc was handled in Eq. (13).
- Confirm the adaptation of [3]'s frequency-separated communication derivation to the shared-band model.
- Obtain original Fig. 3 numerical data or simulation code. No digitized points are used, and visual resemblance is not exact agreement.
