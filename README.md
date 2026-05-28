# NuRNA

> Nussinov Algorithm for RNA Secondary Structure Prediction across Families

<p align="center">
<a href="/src/main.ipynb">Notebook</a>
 |
<a href="/doc">Project Report</a>
</p>

---

## Table of Contents

* [NuRNA](#nurna)
* [About](#about)
* [Project Structure](#project-structure)
* [Requirements](#requirements)
* [How to Run](#how-to-run)
* [Task Assignment](#task-assignment)

---

## NuRNA

_IF3211 Domain-Specific Computation (K01)_
<table>
      <tr align="left">
        <td><b>NIM</b></td>
        <td><b>Name</b></td>
        <td align="center"><b>GitHub</b></td>
      </tr>
      <tr align="left">
        <td>13523004</td>
        <td>Razi Rachman Widyadhana</td>
        <td align="center">
          <div style="margin-right: 20px;">
          <a href="https://github.com/zirachw"><img src="https://github.com/zirachw.png" width="48px;" alt=""/> <br/> <sub><b> @zirachw </b></sub></a><br/>
          </div>
        </td>
      </tr>
      <tr align="left">
        <td>13523002</td>
        <td>Refki Alfarizi</td>
        <td align="center">
          <div style="margin-right: 20px;">
          <a href="https://github.com/l0stplains"><img src="https://github.com/l0stplains.png" width="48px;" alt=""/> <br/> <sub><b> @l0stplains </b></sub></a><br/>
          </div>
        </td>
      </tr>
      <tr align="left">
        <td>13523061</td>
        <td>Darrel Adinarya Sunanda</td>
        <td align="center">
          <div style="margin-right: 20px;">
          <a href="https://github.com/Darsua"><img src="https://github.com/Darsua.png" width="48px;" alt=""/> <br/> <sub><b> @Darsua </b></sub></a><br/>
          </div>
        </td>
      </tr>
      <tr align="left">
        <td>13523009</td>
        <td>Muhammad Hazim Ramadhan Prajoda</td>
        <td align="center">
          <div style="margin-right: 20px;">
          <a href="https://github.com/SayyakuHajime"><img src="https://github.com/SayyakuHajime.png" width="48px;" alt=""/> <br/> <sub><b> @SayyakuHajime </b></sub></a><br/>
          </div>
        </td>
      </tr>
</table>

---

## About

<div align="right">(<a href="#table-of-contents">back to top</a>)</div>

<p align="justify">
This repository is created as part of <b>IF3211 Domain-Specific Computation 2025/2026</b> at <a href="https://itb.ac.id" target="_blank">Institut Teknologi Bandung</a>, Topic 4: DNA &amp; RNA Sequences. The project evaluates the Nussinov dynamic programming algorithm for RNA secondary structure prediction using the ArchiveII dataset, with per-family accuracy analysis using sensitivity, PPV, and F1 metrics.
</p>

<p align="justify">
The primary dataset is <a href="https://huggingface.co/datasets/multimolecule/archiveii"><b>ArchiveII</b></a>, a collection of experimentally verified RNA secondary structures spanning 10 RNA families. All data loading and preprocessing uses the <code>datasets</code> library (HuggingFace) and NumPy exclusively, without external bioinformatics library dependencies for core computation.
</p>

<p align="justify">
The Nussinov algorithm implementation covers the O(n²) DP matrix construction and recursive traceback to reconstruct base pairs in dot-bracket notation, with O(n³) time complexity. Evaluation follows the per-family benchmarking methodology of Mathews (2019), producing accuracy tables per RNA family alongside analysis of GC content and sequence length effects on F1. Visualizations include arc diagrams, DP matrix heatmaps, and biological 2D layouts via <i>forgi</i>.
</p>

---

## Project Structure

<div align="right">(<a href="#table-of-contents">back to top</a>)</div>

```
NuRNA/
  dataset/
  doc/
  src/
    NuRNA/
      __init__.py
      constants.py
      data.py
      nussinov.py
      evaluation.py
      visualization.py
      analysis.py
    .streamlit/
      config.toml
    app.py
    examples.json
    results.json
    main.ipynb
  pyproject.toml
```

---

## Requirements

<div align="right">(<a href="#table-of-contents">back to top</a>)</div>

* Python 3.13+
* [uv](https://docs.astral.sh/uv/) (recommended package manager)

### Core dependencies

* `datasets` for loading ArchiveII from HuggingFace
* `numpy` for DP matrix computation and evaluation metrics
* `pandas` for per-family result tabulation
* `matplotlib` for accuracy plots and arc diagrams
* `forgi` for biological 2D RNA structure visualization

---

## How to Run

<div align="right">(<a href="#table-of-contents">back to top</a>)</div>

### Install dependencies

```bash
uv sync
```

This installs all dependencies and the `NuRNA` package as an editable install, so `from NuRNA.nussinov import predict_structure` works in any notebook without path manipulation.

### Open the notebook interactively

```bash
uv run jupyter notebook src/
```

### Execute the notebook non-interactively

```bash
uv run jupyter nbconvert --to notebook --execute src/main.ipynb
```

### Run the web app

```bash
uv run streamlit run src/app.py
```

---

## Task Assignment

<div align="right">(<a href="#table-of-contents">back to top</a>)</div>

| Name | NIM | Task |
| --- | --- | --- |
| Razi Rachman Widyadhana | 13523004 | Project setup, `__init__.py`, `constants.py`, `data.py`, `app.py`. Report: Introduction, Conclusion, final assembly. |
| Refki Alfarizi | 13523002 | `nussinov.py` (core algorithm: `can_pair`, `build_dp_matrix`, `traceback`, `predict_structure`), `main.ipynb`. Report: Methods, Nussinov algorithm and DP formulation. |
| Darrel Adinarya Sunanda | 13523061 | `visualization.py` (`plot_rna_2d`, `plot_arc_diagram`, `plot_dp_matrix`, `plot_comparison`), `examples.json`, `results.json`. Report: Results & Discussion, structure visualization section. |
| Muhammad Hazim Ramadhan Prajoda | 13523009 | `evaluation.py`, `analysis.py`. Report: Results & Discussion, evaluation and per-family analysis section. |

---

<h3 align="center">
NuRNA &copy; 2026 &middot; 13523002 - 13523004 - 13523009 - 13523061
</h3>
