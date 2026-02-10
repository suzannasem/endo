# Investigating Hormonal Influence on Lesion Growth in Endometriosis

This code is supplementary to my Master's thesis at Bryn Mawr College. Building on the model of lesion growth by Claire Miller et. al (2025), I incorporate chemical kinetics of estrogen and estrogen receptors (ERs) to better understand how fluctuations in hormone levels affect immune function.

This is a work in progress and will be updated as I continue my research.

*File guide:*

Miller folder
* Miller_fn.m: Matlab function file for replicating Miller (2025). Contains main ODE systems for endometrial cycle and immune cells.
* Miller_main.m: Matlab main file for replication Miller (2025). Simulates different regimes using ode15s and produces figures as in paper.

Graham folder:
* graham_params.m: Matlab function file that stores parameter values as a struct.
* Graham_fn.m: Matlab function file for replicating Graham (2023). Contains ODE systems for ovulatory cycle, steriod hormones, and pituitary hormones.
* Graham_main.m: Matlab file to run Graham_fn.m using ode15s. Creates a figure of hormones fit to data from McLKeefe.

Citation:

Miller, C., Germano, D. P. J., Chenoweth, A. M., & Holdsworth-Carson, S. (2025). Mathematical modelling of the immune response during endometriosis lesion onset. bioRxiv. [https://doi.org/10.1101/2025.01.20.633967].
