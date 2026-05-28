# BOPBK

**Version 1.0** (Released on 2026-05-28)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20433589.svg)](https://doi.org/10.5281/zenodo.20433589)

[GitHub Repository](https://github.com/baiweiphys/BOPBK/) | [Releases](https://github.com/baiweiphys/BOPBK/releases)

**`BOPBK`** (BO–Product-Bi-Kappa) is a Matlab numerical code developed by Wei Bai, under the supervision of Dr. Huasheng Xie, who also proposed the BO framework. It is designed for analyzing plasma waves and instabilities in both space and laboratory plasmas, specializing in obliquely propagating waves in magnetized, multi-species hot plasmas.

It supports a range of multi-component velocity distributions:
- Anisotropic drift loss-cone product-bi-kappa (Type I and Type II PBK)
- Anisotropic drift loss-cone kappa-Maxwellian (Type I and Type II KM)
- Anisotropic drift loss-cone Maxwellian-kappa (Type I and Type II MK)
- Anisotropic drift loss-cone bi-Maxwellian (BM)
- Hybrid combinations of the above


## Developer & Supervisor

**Wei Bai** (Developer)  
College of Electrical and Power Engineering  
Taiyuan University of Technology  
Taiyuan 030024, China  
Email: baiweiphys@gmail.com  

**Huasheng Xie** (Supervisor, proposed the BO framework)  
ENN Science and Technology Development Co., Ltd., Langfang 065001, China  
& VeloAlpha Technology Co., Ltd., Beijing 100080, China  
Email: huashengxie@gmail.com  


## Benchmark Examples

Four benchmark cases are provided:

1. **Case 1** — R‑, L‑, and P‑mode waves
2. **Case 2** — Whistler instability
3. **Case 3** — Firehose instability
4. **Case 4** — Ion‑cyclotron instability

See `runall.m` to run all cases sequentially.  
Individual 3-D scanning cases `run_case3_firehose3D.m` and `run_case4dot2_EMIC3D.m` are also provided.


## How to Run

### Run All Cases

Run `runall.m` to perform all simulations at once.  
Data are saved to `results/plt_Fig**/pbkData/`.  
After the runs finish, run `plotall.m` to generate figures (exported as `.pdf` and `.fig`).

### Run a Single Case

1. Go to the case directory, e.g. `BO-PBK/examples/case_***/`.
2. Edit `bopbk.in` to set parameters: `J`, `deg` (propagation angles in degrees), `B0` (background field along z), etc.
3. Run `main_bopbk.m`.
4. Results are plotted and saved in the `output` directory.

### Reproduce a Published Figure

Example for Fig. 4 of Lazar et al. (2010):

1. Run `case02.1_whistler_Lazar2010/main_bopbk.m`, then run `selected-plot/bo_wpdat.m`.
2. Copy the generated `*.dat` files to `../../plt_Fig4/BOPBK_data`.
3. Repeat for `Whistler_Lazar10Fig2_pbk2/` and `Whistler_Lazar10Fig2_pbk6/`, each time copying the `*.dat` files to the same destination.
4. Run `plt_Fig4/plt_benchmark_lazar2010Fig2.m` to produce the figure.

### Note

The MATLAB scripts `run_case3_firehose3D.m` and `run_case4dot2_EMIC3D.m` perform 3-D parameter scans and are time-consuming to run. For efficient computation, the Julia version (BOPBK.jl) is recommended.


## References

Wei Bai, and Huasheng Xie. “BO-PBK: A comprehensive solver for dispersion relations of obliquely propagating waves in magnetized multi-species plasma with anisotropic loss-cone drift product-bi-kappa distribution.” arXiv preprint arXiv:2512.06901, 2025. https://arxiv.org/abs/2512.06901.


## License

This project is licensed under the BSD 3-Clause License — see the [LICENSE](LICENSE) file for details.

---

Created by Wei Bai on 2026-05-25  
