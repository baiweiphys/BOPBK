# BO-PBK

**Version 1.0 beta** | [GitHub Repository](https://github.com/baiweiphys/BOPBK/)


## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.


## Authors & Affiliations

**Wei Bai**  
College of Electrical and Power Engineering  
Taiyuan University of Technology  
Taiyuan 030024, China  
Email: baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn

**Huasheng Xie**  
ENN Science and Technology Development Co., Ltd.  
Hebei Key Laboratory of Compact Fusion  
Langfang 065001, China  
Email: huashengxie@gmail.com, xiehuasheng@enn.cn

## Release Date
2025-12-08

## Program Description


**BO-PBK** (BO-Product-Bi-Kappa), developed by **Wei Bai** and **Huasheng Xie**, is a numerical code for analyzing waves and instabilities in obliquely propagating, magnetized multi-species hot plasmas in space and laboratory settings.


The code supports a wide range of multi-component velocity distributions, such as:
- Anisotropic drift loss-cone product-bi-kappa (PBK)
- Anisotropic drift loss-cone kappa-Maxwellian (KM)
- Anisotropic drift loss-cone bi-Maxwellian (BM)
- Hybrid combinations of these distributions


## Benchmark Examples

Four benchmark cases are provided as representative examples:
1. R-, L-, and P-mode waves
2. Whistler instability
3. Firehose instability
4. EMIC waves


## Development Status

### Current Release
- **MATLAB**: Available (R2024a-compatible)

### In Development
- **Julia**: In progress (porting & migration done) 


## File Descriptions


### Example Cases
- **`examples/Case_XXX/`** - Individual test case directories containing specific configuration files.
  - **`main_bopbk.m`** - Primary driver routine that initiates the simulation with user-defined parameters.
  - **`output/`** - Directory for storing simulation outputs, including dispersion relations and growth rates.
  - (Optional module) **`selected_plot/`** - Contains enhanced visualization routines for 2D/3D plots (based on Xie 2019, CPC).

### Core Source Code
- **`src/`** - Directory containing the core numerical routines and functions.
  - **`solver_pbk.m`** - Numerical solver for plasmas described by PBK distributions.
  - **`solver_maxwell.m`** - Solver for BM plasma models.
  - **`solver_mixed.m`** - Solver for hybrid plasma models combining PBK, KM, and BM distributions.
  - **`func_Jpole.m`** - Implements the J-pole expansion technique for efficient evaluation of plasma dispersion functions.
  - **`getIndexOfBlkMatrix_*.m`** - Functions for extracting indices of sub-block matrices in the global dispersion matrix of PBK, BM, and their mixed.
  - **`maxwell_*.m`** - Coefficient calculation routines for the oblique BM plasma model.
  - **`pbk_*.m`** - Coefficient calculation routines for the oblique PBK plasma model.
  - **`M_*.m`** - Functions that construct equivalent matrices for PBK, BM, and mixed distribution models.
  - **`getLen_*.m`** - Utility functions for determining matrix dimensions and sizes.


## How to run the code:

1. **Run a single test case**:
   - Navigate to the specific case directory: `BO-PBK/examples/case_***/`
   - Configure the input parameters in the input file `bopbk.in`
   - Set parameters such as `J`, `deg` (propagation angles in degrees), and `B0` (background magnetic field in z-direction)
   - Run the main script `main_bopbk.m`
   - Results will be displayed as plots and saved in the `output` directory
2. **Visualize and save data**:
   - Navigate to `selected_plot/`
   - Configure `bo_wpdat.m` and run it to generate figures and save data as `*.dat` files.

3. Show the benchmark cases
   - **Example 1:**
        - 1. Run `case02.1_whistler/main_bopbk`
  Then run `selected-plot/bo_wpdat`
        - 2. Copy `*.dat` files to `../../plt_Fig4/BOPBK_data`
        - 3. Repeat the above steps for `Whistler_Lazar11Fig2_pbk2/` and `Whistler_Lazar11Fig2_pbk6/`. After running each, copy the corresponding `./selected-plot/*.dat` files to `plt_Fig4/BOPBK_data`
        - 4. Run `/plt_Fig4/plt_benchmark_lazar2011Fig2.m`
  
   - **Example 2:**
        - 1. Run `case03.1_firehose/Astfalk2017Fig1_OFHI3d_bm/main_bopbk.m`
        - 2. Copy `/output/*.dat` files to `../../plt_Fig8/pbk_data`
        - 3. Repeat the above steps for `Astfalk2017Fig1_OFHI3d_pbk210main_bopbk.m/`. After running it, copy the corresponding `output/*.dat` files to `plt_Fig8/pbk_data`
        - 4. Run `/plt_Fig8/plt_Astfalk2017fig1_contour2.m`


## References
Wei Bai, and Huasheng Xie. “BO-PBK: A comprehensive solver for dispersion relations of obliquely propagating waves in magnetized multi-species plasma with anisotropic loss-cone drift product-bi-kappa distribution.” arXiv preprint arXiv:2512.06901, 2025, https://arxiv.org/abs/2512.06901.




If you encounter any issues, kindly get in touch with us.

Wei Bai, baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn, TYUT

Huasheng Xie, huashengxie@gmail.com, xiehuasheng@enn.cn, ENN


Created by Wei Bai on 2024-07-24  
Updated by Wei Bai on 2025-12-11


