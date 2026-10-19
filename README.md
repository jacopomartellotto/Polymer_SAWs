<div align="center">

# 🧵 Two-dimensional self-avoiding walks

**Direct rejection and the Pivot algorithm**

*Two sampling strategies. One polymer model.*

🔲 2D square lattice · 🎲 Monte Carlo · ⚛️ Polymer physics

[📊 Results](#results) &nbsp;·&nbsp; [🎓 Presentation](#presentation) &nbsp;·&nbsp; [🚀 Run the project](#run-the-project) &nbsp;·&nbsp; [📚 References](#references)

**Jacopo Martellotto** · Bachelor's Degree in Physics, University of Pisa  
Supervisor: **Prof. Claudio Bonati**

</div>

---

A self-avoiding walk never visits the same lattice site twice. This simple
constraint captures excluded volume in a polymer: two monomers cannot occupy
the same position. This thesis compares two ways of sampling these walks and
examines how their size grows with chain length.

| 🧭 Simple sampling | 🔄 Pivot algorithm |
|---|---|
| Grow independent walks; reject a trial at the first repeated site. | Transform part of an existing walk; reject updates that create overlaps. |
| Direct access to short chains. | Access to much longer chains through a Markov chain. |

<a id="results"></a>

## 📊 Results

<p align="center">
  <a href="slides/figures/01_estensione.png"><img src="slides/figures/01_estensione.png" width="49%" alt="Chain extension and power-law fits from the presentation"></a>
  <a href="slides/figures/02_success_probability.png"><img src="slides/figures/02_success_probability.png" width="49%" alt="Sampling success probabilities from the presentation"></a>
</p>

The figures and values are those shown in the presentation. Chain extension is
fitted with $\langle R_e^2\rangle=AN^{2\nu}$; the theoretical exponent in two
dimensions is $\nu=3/4$.

| Method | Fitted lengths | $\nu$ | $\chi^2_{\mathrm{red}}$ (ndof) |
|---|---|---|---|
| Simple sampling | 6–30 | $0.724\pm0.002$ | 1.53 (12) |
| Pivot | 3,000–20,000 | $0.74996\pm0.00003$ | 1.48 (9) |

Errors represent one standard deviation. Simple sampling gives an effective
short-chain exponent; Pivot agrees with the asymptotic theoretical value.
Success means completing an entire walk for Simple sampling and accepting one
proposed transformation for Pivot.

**⏱️ Run time:** approximately 5 days and 16 hours.

<details>
<summary><strong>🔎 Explore the stability check and the supplementary fit</strong></summary>

**Exponent stability.** Repeat the power-law fit while increasing the shortest
included chain length. The shaded region marks the selected constant-fit range.

![Exponent stability versus the minimum fitted length](slides/figures/09_nu_vs_lmin.png)

**Finite-size corrections.** The backup slides also show the supplementary model
$AN^{2\nu}(1+b/N+c/N^{3/2})$: $\nu=0.75005\pm0.00003$,
$\chi^2_{\mathrm{red}}=1.02$ (17 ndof).

![Supplementary Pivot fit with finite-size corrections](slides/figures/backup_pivot_corrections.png)

</details>

<a id="presentation"></a>

## 🎓 Presentation

<p align="center">
  <a href="slides/tesi_saw.pdf"><img src="slides/preview.png" width="580" alt="Click the cover to open the thesis slides"></a>
</p>

<p align="center">
  <a href="slides/tesi_saw.pdf"><strong>▶ View the slides</strong></a> &nbsp;·&nbsp;
  <a href="slides/tesi_saw.tex">Beamer source</a> &nbsp;·&nbsp;
  <a href="long_abstract.pdf">Long abstract</a>
</p>

<a id="run-the-project"></a>

## 🚀 Run the project

From PowerShell in the project folder, create the Python environment and install
the pinned dependencies:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

**Reanalyse saved measurements and regenerate outputs:**

```powershell
.\.venv\Scripts\python.exe -B scripts/render_outputs.py
```

**Run the full simulation pipeline** — also requires `g++` with C++17 support:

```powershell
.\run_all.ps1 -Python ".\.venv\Scripts\python.exe"
```

The pipeline calibrates Pivot thermalization, runs production, and generates
fits, tables and figures. Completed compatible replicas are reused when
resuming. Seeds and simulation settings are recorded in the protocol.
Power-law fits use `scipy.optimize.curve_fit` with `absolute_sigma=True`.

To rebuild the slide PDF with a LaTeX installation, run `.\slides\build.ps1`.

## 🗂️ Find your way around

| Location | What you will find |
|---|---|
| [config/](config/) | Lengths, seeds, replica counts and thermalization settings |
| [grid/](grid/) · [src/](src/) | Simple-sampling and Pivot simulation code |
| [scripts/](scripts/) | Production, analysis and plotting |
| [outputs/figures/](outputs/figures/) | Generated PNG plots and comparison table |
| [outputs/data/](outputs/data/) · [outputs/analysis/](outputs/analysis/) | Measurements, fit parameters and uncertainties |
| [outputs/diagnostics/](outputs/diagnostics/) | Thermalization traces and production cutoffs |
| [slides/](slides/) | Presentation, source and the figures displayed in this README |

<a id="references"></a>

## 📚 References

- P.-G. de Gennes, *Scaling Concepts in Polymer Physics*. Cornell University Press (1979).
- K. Binder and D. W. Heermann, *Monte Carlo Simulation in Statistical Physics*, 4th ed. Springer (2002). [DOI](https://doi.org/10.1007/978-3-662-04685-2).
- N. Madras and A. D. Sokal, *J. Stat. Phys.* **50**, 109–186 (1988). [DOI](https://doi.org/10.1007/BF01022990).
- S. Caracciolo et al., *J. Stat. Phys.* **120**, 1037–1100 (2005). [DOI](https://doi.org/10.1007/s10955-005-7004-3).
- R. Bauerschmidt et al., “Lectures on self-avoiding walks”, *Clay Math. Proc.* **15** (2012). [Lecture notes](https://www.ihes.fr/~duminil/publi/saw_lecture_notes.pdf).
