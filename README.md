# **Periodic component analysis ($\pi$ CA)**

Extraction of periodic components from multidimensional recordings

The methods are described in [this article](https://link.springer.com/chapter/10.1007/978-3-319-93764-9_48)

## Description

This repository contains the codes written to extract periodic components from multidimensional recordings using spatial filtering. The period of interest $T = 1/F$ should be known. 

All codes are written in Matlab. 

If we denote the filtered signal by $s(t) = w^t \mathbf{x}(t)$, the filters are obtained as:

Method 1: 

![equation](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bw%7D_*=%5Carg%5Cmin_%7B%5Cmathbf%7Bw%7D%7D%5Cfrac%7B%5Csum_t(s(t&plus;T)-s(t))%5E2%7D%7B%5Csum_t(s(t))%5E2%7D)

Method 2: 

![equation](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bw%7D_*=%5Carg%5Cmax_%7B%5Cmathbf%7Bw%7D%7D%5Cfrac%7B%5Csum_t%20s(t&plus;T)%5Ccdot%20s(t)%7D%7B%5Csum_t(s(t))%5E2%7D)

Method 3 (canonical correlation analysis): 

![equation](https://latex.codecogs.com/svg.image?(%5Cmathbf%7Bw%7D_*,%5Cmathbf%7Bw%7D_y_*)=%5Carg%5Cmax_%7B%5Cmathbf%7Bw%7D,%5Cmathbf%7Bw%7D_y%7D%5Cfrac%7B%5Cmathbf%7Bw%7D%5ET%20C_%7Bx;y%7D%5Cmathbf%7Bw%7D_y%7D%7B%5Csqrt%7B%5Cmathbf%7Bw%7D%5ET%20C_%7Bx%7D(0)%5Cmathbf%7Bw%7D%5Ccdot%5Cmathbf%7Bw%7D_y%5ET%20C_%7By%7D(0)%5Cmathbf%7Bw%7D_y%7D%7D)

with $y(t) = (sin(2\pi F t), cos(2\pi F t), sin(2\pi 2F t), \ldots, sin(2\pi N_h F t), cos(2\pi N_h F t))^T$ ($N_h$ being a parameter: the number of harmonic of the reference periodic signal) and $C_{x;y} =\mathbb{E}_t\[x(t)y(t)^T\]$. 

Method 4 (Spectral Contrast maximization): 

![equation](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bw%7D_*=%5Carg%5Cmax_%7B%5Cmathbf%7Bw%7D%7D%5Cfrac%7B%5Cmathbb%7BE%7D_%7Bf%5Cin%5Cnu%7D%5C%7B%7CS(f)%7C%5E2%5C%7D%7D%7B%5Cmathbb%7BE%7D_%7Bf%5Cin%5Cmu%7D%5C%7B%7CS(f)%7C%5E2%5C%7D%7D)

with the Fourier transform of the filtered signal at frequency $f$: $S(f) := \mathcal{F}_f\[s(t)\] = \mathbf{w}^T \mathcal{F}_f\[\mathbf{x}(t)\] = \mathbf{w}^T X(f)$


## Running

An example can be run using test_piCA.m

The methods are implemented in piCA_compute.m

![Image](https://github.com/user-attachments/assets/15159a7e-aae9-4c8d-b42a-bd3c41326977)

## Contact

You can contact me at dounia **dot** mulders **at** uclouvain.be for any question. :-)


