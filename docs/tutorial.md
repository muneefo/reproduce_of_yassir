# Reading the physical system

Source: the complete five-page arXiv:2604.06069v1 PDF supplied by the user (7 April 2026). Communication derivations: arXiv:2502.13348v2 (4 December 2025), §IV and Appendices F–G, corresponding to reference [3]. The latter is a preprint; equality to the final journal version has not been established. Neither PDF is redistributed in this repository.

## What transmits, what receives, and where ISAC occurs

BS locations are a planar homogeneous Poisson point process (PPP). A BS makes eight directional beams. Each beam points toward one of its own sensing targets and serves a scheduled mobile user in that direction. ISAC occurs **within each beam and slot**, using a common known signal structure and shared energy. The user is assumed plentiful enough that every beam has someone to serve. This is the original terrestrial model; there is no satellite, spherical geometry, altitude, or orbital motion here.

At the start of a slot the BS sends a sensing pulse of duration Ts=1/W. It allocates γE joules to that pulse. Over the remaining Tt−Ts it sends communication data and receives delayed echoes simultaneously. Thus Ps=γE/Ts and Pc=(1−γ)E/(Tt−Ts). This is not a long sensing interval followed by an independent communication-only interval. Full-duplex reception during downlink transmission needs cancellation of the BS's own transmit signal. The residual contribution is Pc ζSI.

For one target, its nearest BS is the sensing receiver. Its **monostatic** echo travels BS → target → same BS. A neighboring BS's **bistatic** echo travels neighboring BS → target → sensing BS. The neighboring BS continues its existing service; it does not steer an extra beam for this target. It contributes only if one of its existing beams happens to illuminate the target over a LoS link. Membership in the nearest-Nc cluster alone does not guarantee contribution.

The sensing BS adds these echo **powers**, not their complex amplitudes. This non-coherent combination requires neither carrier-phase alignment nor coherent joint beamforming. It nevertheless requires a common known transmitted structure, synchronization, and timing coordination sufficient to identify/exploit the echoes. The paper assumes that connectivity supports this coordination and that the receiver can exploit the echoes. It does not specify a practical waveform-separation algorithm, its estimation errors, processing cost, or timing tolerance. No remote sensing observations are sent to a fusion center; this sensing BS processes its received echoes locally.

There are two network interference paths. **Direct interference** travels from another BS straight to the sensing receiver, if transmit and receive beams align. It can be LoS or NLoS and undergoes one-hop propagation and Nakagami fading. **Target-reflected interference** travels from an out-of-cluster BS via this target to the sensing receiver; it is an opportunistically illuminated LoS echo that has not been classified as useful. The same BS can contribute to both terms. Their geometries and beam directions are therefore related. Cluster BSs are removed only from reflected interference, not automatically from direct interference.

Blockage of a nonserving link has LoS probability exp(−ηb r). The serving target and user links are forced LoS. Echo fluctuations are independent exponential radar-cross-section (RCS) draws, using the Swerling I assumption. We do not multiply useful echoes by additional Nakagami gains: the printed radar powers already specify exponential RCS variation, while Nakagami fading appears in the one-hop direct and communication links.

## Parameter and notation table

Values below come from §IV unless identified as derived or an implementation choice. SI values are used in `config.py`. λ without a subscript in the radar equation is wavelength, not PPP density.

| Symbol | Meaning | Internal value/unit | Source / implementation |
|---|---|---|---|
| ΦBS, λBS | BS PPP and intensity | 70/10⁶ m⁻² | §II, Eq. (15); `Config.density_m2` |
| R1,…,RNc | Ordered BS distances to target | m, random | Eqs. (12),(15); `ordered_radii`, `distance_qmc` |
| ri | Interferer distance to receiver BS | m, random; distinct from Ri | Eqs. (9),(10); `spatial_marks` |
| RI | Printed direct-interference guard | m, random | Eq. (11); `DirectLT` |
| Nc | Sensing cluster size | 1,2,4,6; also 8 for density | Fig. 3; 1 is monostatic |
| Nb | Simultaneous beams per BS | 8 | §IV |
| W | Shared bandwidth | 200×10⁶ Hz | §IV |
| f, λ | Carrier / wavelength | 28×10⁹ Hz; c/f≈0.010707 m | §IV / radar equation |
| θb | 3-dB full beamwidth | 18π/180 rad | Eqs. (1),(4) |
| Gm | Maximum **power** antenna gain | 10^(20/10)=100 | 20 dBi; Eqs. (1),(8)–(10) |
| q | Cosine spread | π/θb=10, derived | Eq. (1); half-power at ±θb/2 |
| αL,αN | Path exponents | 2,4; dimensionless | Eqs. (8)–(10) |
| CL,CN | Reference-distance power intercepts | 10^(−61.4/10),10^(−72/10) | Eqs. (9),(10),(21), referenced to 1 m |
| mL,mN | Nakagami power shape | 3,2 | Mean-one Gamma(m,1/m) in simulation |
| σav,mono,σav,bi | Mean radar cross sections | 1,0.7 m² | Unnumbered radar equations after (4) |
| ηb | Blockage coefficient | 0.0149 m⁻¹ | Eq. (2); units implied by exponent |
| γ | Sensing energy fraction | 0.7 | §II, §IV; configurable [0,1] |
| Pavg | Slot-averaged power budget | 1 W | §IV; budget scope ambiguous |
| Reff | Effective Voronoi radius | (πλBS)^−1/2≈67.43 m | §IV |
| Rmax | Nominal unambiguous range | 3Reff≈202.30 m | §IV; **not** an imposed distance cutoff in equations |
| c | Propagation speed | 299792458 m/s | Physical constant implementation choice; paper leaves c symbolic |
| Ts,Tt | Pulse / slot durations | 5 ns / ≈1.3496 μs | §II |
| E | Slot energy | Pavg Tt≈1.3496 μJ | Derived; per-beam default |
| Ps,Pc | Pulse / communication power | ≈188.945 W / ≈0.301116 W | Derived |
| T,KB | Temperature / Boltzmann constant | 300 K; 1.380649×10⁻²³ J/K | §IV / physical constant |
| KB T W | Noise power | 8.283894×10⁻¹³ W | Eq. (7); no noise figure specified |
| ζSI | Residual power fraction | 10⁻¹² | Eq. (7); 120 dB power suppression |
| τ,t,b | SINR threshold, reliability, moment order | dimensionless; b may be complex | Eqs. (14),(18)–(20) |
| keff,θeff | Gamma matched shape / scale | dimensionless / W | Eqs. (5),(6) |
| η | Alzer coefficient | Γ(1+keff)^(-1/keff) | Eq. (14); distinct from ηb |
| sz | Laplace argument | W⁻¹ | zητ/θeff in Eq. (14) |
| area, networks | Paper simulation setting | 25×10⁶ m²; 200000 | §IV; saved run counts may be smaller |
| seed | Random-number initialization | 20261010 plus documented offsets | Implementation choice, not a paper parameter |

Use 10^(dB/10) for **power ratios**, including SINR, intercepts, and gains. The radar and sensing direct-interference equations require Gm²=10000, not 100 and not 10^(20/20). The radar equation has λ²/(4π)³; do not insert CL a second time. At exponents other than free-space 2, a reference-distance convention would be needed for dimensional consistency. The paper's stated αL=2 avoids that issue for the radar expression.

## Equation-by-equation reading and derivation

**(1): angle reduces gain.** The cosine-squared power pattern peaks at boresight and vanishes outside ±π/q. Solving cos²(qθb/4)=1/2 gives q=π/θb. The support is therefore a full **2θb**, not θb. Sensing analysis explicitly replaces it with a top-hat gain Gm on a full θb support. `beam_offset` implements shared beam orientation; communication uses the cosine, sensing uses the stated top-hat. No sidelobes are added.

**(2): longer paths are more often blocked.** Exponential distance decay pLOS(r)=exp(−ηb r) is used as an independent link mark conditional on geometry. The paper gives marginal probabilities, not a random-building/blocker process, so correlations due to the same physical blocker cannot be reconstructed. Serving links are forced LoS rather than obtained by associating to the nearest LoS BS.

**(3),(4): count only available cooperative echoes.** Write Yn=Bn Xn with Bn an availability indicator. A uniformly rotated non-overlapping beam set covers an angular fraction Nb θb/(2π). Multiplying by the LoS probability gives pn. The sum Ptot=X1+ΣYn adds energies/powers without cross-phase terms. `sensing_slot` samples the original echoes and uses the same target mark for useful and interfering classifications.

**Radar powers after (4): two legs create severe attenuation.** Monostatic mean is A σav,mono/R1^(2αL); bistatic mean is A σav,bi/(R1^αL Rn^αL), with A=Ps Gm² λ²/(4π)³. The same receiver-target return leg R1 occurs for all echoes. An exponential RCS remains exponential after multiplication by a deterministic radar factor. There is no extra fading term in these printed powers.

**(5): on/off echoes have more variance than continuously present echoes.** E[Yn]=pn μn; E[Yn²]=pn·2μn². Therefore Var(Yn)=pn(2−pn)μn². Independent useful echo marks add their variances. `echo_moments` implements these exact conditional moments. Independent equally spaced rotations and independent target-link LoS marks give these marginals; their dependence on the direct-interference beam marks is retained by simulation.

**(6): Gamma fitting preserves only two moments.** For Gamma(k,θ), E[X]=kθ and Var(X)=kθ². Solving yields k=mean²/variance and θ=variance/mean. This is an approximation to the Bernoulli-exponential sum, not its exact law. Monostatic Nc=1 gives k=1 and an exact exponential. We quantify the Gamma CDF error separately from the Alzer error in `validation.json`.

**(7)–(10): SINR is desired echo divided by all disturbances.** The numerator is the original power sum; reflected interference uses the same radar factor with outside BS distances Ri. Direct LoS/NLoS interference uses one-hop intercepts, receive and transmit gains, and Gamma(m,1/m) power fading. Noise and residual SI are deterministic. These definitions, with shared geometry and marks, are implemented in `sensing_slot`. No Gamma approximation is used in simulation.

**(11): average a PPP interference transform over a random guard.** A mean-one Gamma link gives E[exp(−s Ph)]=(1+sP/m)^−m. A PPP product becomes exp(−∫intensity·[1−link-LT]darea) by the PGFL. The two antenna alignment fractions yield ξ=2πλBS Nb(θb/2π)². Integrating r dr from RI and averaging RI with density 2πλBS r exp(−πλBS r²) gives the printed expression. `DirectLT` uses logarithmic radial integration and Gauss–Laguerre integration of πλBS RI²~Exp(1). This printed guard is an approximation; see the audit.

**(12): use squared-area arrivals for PPP neighbor distances.** Given R1, πλBS(R_(Nc+1)²−R1²) is a sum of Nc independent Exp(1) arrivals, hence Gamma(Nc,1). A change of variables produces the printed density. It is correct given only R1. Given the entire vector R, however, the next gap is exponential after RNc², not after R1² with shape Nc. `reflected_lt` deliberately follows the printed choice.

**(13): a marked exponential echo has a rational LT.** A target-interference mark contributes 1/(1+sAσav,bi/(R1^αL r^αL)), with illumination probability folded into λBS Nb θb exp(−ηb r). Applying PGFL and averaging the printed random radius gives (13). For αL=2, `reflected_integral` evaluates the radial integral through the complex exponential integral E1. The factor r/(r²+a) splits into two simple fractions with poles ±i√a, giving a·Re[exp(−iηb√a)E1(ηb(lower−i√a))]. Independent quadrature validates this reduction. The omitted nearest-point contribution is documented, not repaired in the printed curves.

**(14): replace a Gamma survival function by an exponential expansion.** For x=τItot/θeff, the Alzer survival approximation is 1−(1−exp(−ηx))^keff. A generalized binomial expansion expresses it as Σz≥1 cz exp(−zηx), cz=(−1)^(z+1) binom(keff,z). Averaging each exponential yields the product of interference LTs and deterministic noise/SI terms. `conditional_coverage` uses cz+1=−cz(keff−z)/(z+1); it never rounds keff to an integer. After z≥ceil(keff), all remaining coefficients have one sign, and the omitted coefficient mass is |cz|(z−keff)/keff. Multiplication by the current decreasing LT gives a conservative absolute probability-error bound. Failure to attain the chosen tolerance raises an error. The paper's **independence** step is distinct from this numerical series approximation.

**(15),(16): average over where the target's neighboring BSs lie.** A PPP has no points inside RNc except the ordered Nc points; combining its void probability with Nc radial densities yields (2πλBS)^Nc(∏rn)exp(−πλBS rNc²). Simulation and scrambled Sobol quadrature generate equivalent ordered radii from independent exponential squared-area gaps, avoiding a difficult ordered multidimensional integral. `distance_qmc` is spatial integration, not a change of model. Its finite sample error is assessed by increasing its size.

**(17): integrate the SINR tail to get an information rate.** For nonnegative X, E[log(1+X)]=∫0∞P(X>exp(v)−1)dv. Multiplying by Nb counts beams per BS. `sensing_rate` integrates analytical coverage; the Monte Carlo estimator is Nb·mean(log1p(SINR)). Both use natural logarithms, hence nats/s/Hz. This is a sensing information metric, not an actual decoded data throughput. Summing it with communication rates is the paper's comparison metric, not a proof of additive operational capacity.

**(18),(19): average success and reliability dispersion answer different questions.** PS(R) averages over everything not specified by ordered cluster distances. Its meta-distribution is the distribution of that conditional probability across distance vectors. Mb=E[PS(R)^b], so M1 must equal average coverage for the same model and conditioning. For example, targets with success probabilities 0.99 and 0.61 in equal proportions have average success 0.80, but only 50% exceed reliability 0.90. A population with every target at 0.80 has the same average and **zero** exceeding 0.90. `moments` and nested distance simulation use the paper's conditioning; the full-marked-geometry comparison is deliberately labeled differently.

**(20): Fourier-invert log reliability.** Mju=E[exp(ju log PS)] is a characteristic function of log PS. Gil–Pelaez inversion gives its CCDF at log t. `gil_pelaez` handles the u→0 limit analytically, then integrates a finite frequency range. It returns raw values even if truncation causes ringing. The empirical CCDF of the same PS quadrature samples supplies an independent numerical inversion comparison, but not an independent model validation.

**(21),(22): communication has one useful hop and an omnidirectional user receiver.** Desired power is Pc Gm CL R0^−αL times mean-one Nakagami power fading, with no residual SI at the user. Nearest-BS association excludes interferers closer than R0. Target reflections are not added separately to user SINR; [3] treats their multipath effects as part of channel fading. [3, Eqs. (24)–(26)] averages the cosine pattern over uniform interfering beam directions and uses the finite mL=3 Alzer sum with kL=mLΓ(1+mL)^−1/mL. `communication.py` implements that derivation and a clearly identified Nb-fold shared-band adaptation. The prefactor Nb(Tt−Ts)/Tt in (22) counts user beams and the communication fraction of a slot.
