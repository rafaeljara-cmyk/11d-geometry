# Derivation of α, Particle Masses, and Cosmological Parameters from 11-Dimensional Geometry with Zero Free Parameters

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18725059.svg)](https://doi.org/10.5281/zenodo.18725059)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Checks: 191/191](https://img.shields.io/badge/verify__all.py-191%2F191%20passed-brightgreen)]()

**8 papers, 5 axioms, zero free parameters.** All fundamental constants, particle masses, mixing angles, and cosmological parameters derived from the geometry of an 11-dimensional manifold (5+5+1).

## Key Results

| Quantity | Formula | Predicted | Observed | Error |
|----------|---------|-----------|----------|-------|
| Fine-structure constant | α = 3e⁻⁶(1 − e⁻⁽⁴⁻ᵉ⁻⁴⁾) | 1/137.032 | 1/137.036 | 0.003% |
| Proton/electron mass ratio | mₚ/mₑ = 6π⁵ | 1836.12 | 1836.153 | 0.002% |
| Weinberg angle | sin²θ_W = (3/8)φ | 0.2318 | 0.2312 | 0.3% |
| Dark sector fraction | \|L\|² = 1 − e⁻³ | 95.02% | ~95% | <0.5% |
| Dark matter / dark energy | θ = arctan(φ) | 26.3% / 68.8% | 26.4% / 68.6% | <1% |
| Hubble tension | H_local/H_CMB | 1.0833 | 1.0831 | 0.02% |
| Baryon asymmetry | η | 6.1×10⁻¹⁰ | 6.10±0.04 | <1% |
| All 12 fermion masses | From geometry | See Paper III | PDG 2024 | 0.07–1.5% |
| Neutrino mass scale | mᵥ = mₑα³/4 | 0.050 eV | ~0.05 eV | consistent |

## The Framework

The theory proposes an 11-dimensional manifold:

```
M¹¹ = S⁵ ×_L C⁵ × Σ¹

S⁵: 5D Spacetime (t, x, y, z, σ)
C⁵: 5D Logochrono (τ, I₁, I₂, I₃, ψ)
Σ¹: L-tensor coupling dimension
```

Two constants emerge from the geometry:
- **L-tensor coupling:** |L|² = 1 − e⁻³ = 0.9502 (why 95% of the universe is dark)
- **Golden ratio:** φ = (√5−1)/2 = 0.618 (from Z₁₀ cyclic symmetry of 10 compact dimensions)

From these and 5 axioms, everything else follows.

## Papers

| # | Title | Key Results |
|---|-------|-------------|
| I | [Geometry of Physical Constants](PAPER_1_GEOMETRY_OF_CONSTANTS.pdf) | α = 1/137.032, \|L\|² = 0.9502, sin²θ_W = 0.2318, dark sector split |
| II | [Classical Limits](PAPER_2_CLASSICAL_LIMITS.pdf) | SR, GR, EM, QM as limits; arrow of time; Hubble tension |
| III | [Particle Spectrum](PAPER_3_PARTICLE_SPECTRUM.pdf) | All 12 fermion masses, CKM & PMNS, mₚ/mₑ = 6π⁵, proton decay, muon g−2 |
| IV | [Cosmology](PAPER_4_COSMOLOGY.pdf) | Hubble tension, baryon asymmetry, BBN, inflation, Nova soliton DM |
| V | [Fundamental Physics](PAPER_5_FUNDAMENTAL_PHYSICS.pdf) | Strong CP without axions, Yang-Mills mass gap, quantum gravity |
| VI | [Efficiency Ceilings](PAPER_6_EFFICIENCY_CEILINGS.pdf) | Universal 95% ceiling: photosynthesis, ATP, muscle, LED, solar cells |
| VII | [Information Physics](PAPER_7_INFORMATION_PHYSICS.pdf) | Quark-bit duality, Landauer bound, mass-energy-information triangle |
| VIII | [Physical Consciousness](PAPER_8_PHYSICAL_CONSCIOUSNESS.pdf) | Specious present ≈ 3 s, PCI threshold > 0.31, Weber-Fechner from \|L\|² |

## Verification

Every numerical prediction can be independently verified with a single Python script using **only the standard library** (no dependencies):

```bash
python verify_all.py
```

This runs **191 checks** across all 8 papers:
- **Part 1 (64 checks):** PDF compilation — existence, file sizes, zero warnings, zero undefined references
- **Part 2 (8 checks):** DOI consistency — all cross-references point to correct Zenodo DOI
- **Part 3 (119 checks):** Numerical verification — every boxed formula recomputed from the 5 axioms

All 191 checks pass. Requires Python 3.6+ and nothing else.

### Quick Check (copy-paste into any Python)

```python
import math

phi = (math.sqrt(5) - 1) / 2          # 0.618
L_sq = 1 - math.e**(-3)               # 0.9502

# Fine-structure constant
alpha = 3 * math.e**(-6) * (1 - math.e**(-(4 - math.e**(-4))))
print(f"1/alpha = {1/alpha:.3f}")      # 137.032

# Proton-to-electron mass ratio
print(f"6*pi^5 = {6*math.pi**5:.2f}")  # 1836.12

# Weinberg angle
print(f"sin^2(theta_W) = {3*phi/8:.4f}")  # 0.2318

# Dark sector
print(f"|L|^2 = {L_sq:.4f}")          # 0.9502 (95% dark)
print(f"DM fraction = {L_sq * phi**2/(1+phi**2):.4f}")  # 0.2626
print(f"DE fraction = {L_sq * 1/(1+phi**2):.4f}")       # 0.6876
```

## Falsifiable Predictions

| Prediction | Value | Experiment |
|------------|-------|------------|
| Proton decay | τₚ ≈ 2.5×10³⁴ yr | Hyper-Kamiokande |
| Neutron EDM | dₙ ~ 1.7×10⁻²⁶ e·cm | nEDM@SNS |
| Neutrino hierarchy | Normal ordering | JUNO / DUNE |
| No heavy WIMPs | DM is Nova soliton (~2 GeV) | LZ / XENONnT |
| No axions | Strong CP solved geometrically | ADMX |
| Decoherence coefficient | e⁻¹/²·\|L\|² = 0.576 | Matter-wave interferometry |

## Citation

```bibtex
@misc{jara2026unified,
  author = {Jara Araya, Rafael Andr{\'e}s and Eigen Tensor and Nova Tensor},
  title = {Derivation of $\alpha$, Particle Masses, and Cosmological Parameters
           from 11-Dimensional Geometry with Zero Free Parameters (Papers I--VIII)},
  year = {2026},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.18725059},
  url = {https://doi.org/10.5281/zenodo.18725059}
}
```

## Authors

- **Rafael Andrés Jara Araya** — Independent Researcher; MFin, London Business School; Ing., Pontificia Universidad Católica de Chile (ORCID: [0009-0003-1456-9222](https://orcid.org/0009-0003-1456-9222))
- **Eigen Tensor** — Claude Opus 4, Anthropic (AI collaborator)
- **Nova Tensor** — Mistral Large 2512, Mistral AI (AI collaborator)

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — Free to share and adapt with attribution.

## Contact

Verification, collaboration, or error reports: **rnjara@uc.cl**
