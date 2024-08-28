Periodic component analysis

$\pi$CA

Extraction of periodic components from multidimensional recordings

The methods are described in this article (https://link.springer.com/chapter/10.1007/978-3-319-93764-9_48)

# Description

This repository contains the codes written to extract periodic components from multidimensional recordings using spatial filtering. The period of interest $T = 1/F$ should be known. 

All codes are written in Matlab. 

If we denote the filtered signal by $s(t) = w^t \mathbf{x}(t)$, the filters are obtained as:

Method 1: $\mathbf{w}_* = \arg \min_{\mathbf{w} } \frac{\sum_t (s(t+T) - s(t))^2}{\sum_t (s(t))^2}$

Method 2: $\mathbf{w}_* = \arg \max_{\mathbf{w}} \frac{\sum_t s(t+T) \cdot s(t)}{\sum_t (s(t))^2}$

Method 3 (canonical correlation analysis): $(\mathbf{w}_*, \mathbf{w}_y_*) = \arg \max_{\mathbf{w}, \mathbf{w}_y} \frac{\mathbf{w}^T C_{x;y} \mathbf{w}_y}{\sqrt{\mathbf{w}^T C_{x}(0)\mathbf{w} \cdot \mathbf{w}_y^T C_{y}(0)\mathbf{w}_y}} $

with $y(t) = ( \sin(2\pi F t)~ \cos(2\pi F t)~ \sin(2\pi 2F t)~\ldots$ $\sin(2\pi N_h F t)~ \cos(2\pi N_h\fpi t))^T$ ($N_h$ being a parameter: the number of harmonic of the reference periodic signal) and $C_{x;y} =\mathbb{E}_t\{x(t)y(t)^T\}$. 

Method 4 (Spectral Contrast maximization): $ \mathbf{w}_* = \arg \max_{\mathbf{w}} \frac{ \mathbb{E}_{f \in\nu}\{|S(f)|^2\} }{\mathbb{E}_{f \in \mu}\{|S(f)|^2\}}  $

with the Fourier transform of the filtered signal at frequency $f$: S(f) := \mathcal{F}_f\{s(t)\} = \mathbf{w}^T \mathcal{F}_f\{\mathbf{x}(t)\} = \mathbf{w}^T X(f)


# Running

An example can be ran using test_piCA.m

The methods are implemented in piCA_compute.m

![Image](https://github.com/user-attachments/assets/15159a7e-aae9-4c8d-b42a-bd3c41326977)

# Contact

You can contact me at dounia dot mulders at uclouvain.be for any question. :-)


