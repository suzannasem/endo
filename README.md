# The LESIOn Model for Endometriosis: Local Estradiol Synthesis in Infertility Onset

This code is supplementary to my Master's thesis at Bryn Mawr College. We build on a model of lesion growth and immune dysfunction by Claire Miller et. al (2025) with the reduced Graham model (2023) for folllicular growth and hormonal dynamics. 

We propose a new mathematical model, the LESIOn Model (Local Estradiol Synthesis in Infertility Onset), which incorporates immune and endometrial cell profiles as well as local production of $E_2$ by endometriotic lesions to investigate ovulation-related infertility in endo. Additionally, we integrate estrogen receptor (ER) binding dynamics to analyze the role of $ER_\beta$ in immune system dysfunction, finding a **higher amount of $ER_\beta$ on endometrial cells in disease**, with the **$ER_\beta/ER_\alpha$ ratio correlated to disease severity**. Disease state could be determined by an individual's epigenetic affinity for $ER_\beta$ expression, suggesting that treatments targeting the ESR2 gene may have therapeutic effects in endometriosis. Higher amounts of locally synthesized estradiol were found to **reduce ovarian follicle growth, reduce $P_4$ surges, and lowered the LH/FSH ratio, implicating local estradiol synthesis in endometriosis-associated infertility.**

*File guide:*

Miller folder
* Miller_fn.m: Matlab function file for replicating Miller (2025). Contains main ODE systems for endometrial cycle and immune cells.
* Miller_main.m: Matlab main file for replication Miller (2025). Simulates different regimes using ode15s and produces figures as in paper, with parameter sweep for NK cells ($\mu_K$, $\sigma$)

Graham folder:
* graham_params.m: Matlab function file that stores parameter values as a struct.
* Graham_fn.m: Matlab function file for replicating Graham (2023). Contains ODE systems for ovulatory cycle, steriod hormones, and pituitary hormones.
* Graham_main.m: Matlab file to run Graham_fn.m using ode15s. Creates a figure of hormones fit to data from McLKeefe and figure varying basal $e_0$ production rate.

Constant Model folder:
* Constant_E2_fn: MATLAB function file for constant model (receptors scaled to be dimensionless)
* Constant_E2_params: MATLAB function file that stores parameter values as a struct.
* Constant_E2_main: Main MATLAB for simulating constant model (shows K1 vs K2 heatmaps, ERb/ERa ratios, and immune/endo cell profiles)

Spline Model folder:
* Variable_E2_fn: MATLAB function file for spline model (receptors scaled to be dimensionless)
* Variable_E2_main: Main MATLAB for simulating spline model (shows K1 vs K2 heatmaps, ERb/ERa ratios, and immune/endo cell profiles)

FrankenModel folder:
* Franken_fn: MATLAB function file for frankenmodel (receptors scaled to be dimensionless)
* Franken_params: MATLAB function file that stores parameter values as a struct.
* Franken_main: Main MATLAB for simulating frankenmodel

LESIOn folder:
* LESIOn_fn: MATLAB function file for LESIOn model (receptors scaled to be dimensionless)
* LESIOn_params: MATLAB function file that stores parameter values as a struct.
* LESIOn_main: Main MATLAB for simulating LESIOn model

figs folder
* All relevant figures across models

Thesis_Draft.pdf: Working draft of my thesis! Contains biological background, mathematical explanations, and LESIOn interpretation.

McL_Keefe.dat: Data used for Graham, Constant, and Spline Model (and for comparison in FrankenModel and LESIOn).

Citations for Miller & Graham:

Miller, C., Germano, D. P. J., Chenoweth, A. M., & Holdsworth-Carson, S. (2025). Mathematical modelling of the immune response during endometriosis lesion onset. bioRxiv. [https://doi.org/10.1101/2025.01.20.633967].

Graham, E. J., Elhadad, N., & Albers, D. (2023). Reduced model for female endocrine dynamics: Validation and functional variations. Mathematical biosciences, 358, 108979. https://doi.org/10.1016/j.mbs.2023.108979
