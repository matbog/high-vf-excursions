# Oscillatory covariance and Bochner diagnostics

This exploration has two separate parts: a proposed covariance built from
a complex scale, and a real field generated from an explicitly nonnegative
radial spectrum. They must not be interpreted as the same model.

## Complex-scale candidate

For $C(r)=\exp[-r^2/(2R^2e^{2i\theta})]$, the relevant real-field candidate
is $\operatorname{Re}C(r)$. Its spectral density, up to a positive Fourier
normalization, is

$$
\operatorname{Re}S(k)
\propto e^{-k^2R^2\cos(2\theta)/2}
\cos\!\left(3\theta-\tfrac12 k^2R^2\sin(2\theta)\right),
\qquad 0\leq\theta<\pi/4.
$$

At $\theta=0$ this is Gaussian and nonnegative. For every
$0<\theta<\pi/4$, the cosine takes negative values at some finite
wavenumber: the candidate fails Bochner's criterion on the full continuum.
A plot over $0\leq k\leq\pi/dx$ tests only that band; underflow and finite
cutoffs can hide negative lobes. It is not a covariance-validity phase
boundary. The endpoint $\theta=\pi/4$ is excluded from this integrable
Gaussian-transform calculation.

The complex even function itself also fails the Hermitian covariance
condition for a stationary complex field unless it is real. The diagnostic
therefore explicitly concerns its real part.

## A valid nonnegative-spectrum alternative

The generator instead uses

$$
S(k)=e^{-L_d^2(k-\omega)^2/4}+e^{-L_d^2(k+\omega)^2/4},
\qquad \omega=2\pi/\lambda,\quad k=\|\mathbf k\|.
$$

This radial spectrum is nonnegative. Its normalized **three-dimensional**
continuum covariance is

$$
C(r)=e^{-(r/L_d)^2}
\frac{\cos(\omega r)+b\,\operatorname{sinc}(\omega r)}{1+b},
\qquad b=\frac{L_d^2\omega^2}{2},
$$

where $\operatorname{sinc}(x)=\sin(x)/x$.
The expression follows from the radial Fourier transform
$\int_0^\infty k^2S(k)\operatorname{sinc}(kr)\,dk$, normalized at $r=0$.
The earlier reference $e^{-(r/L_d)^2}\operatorname{sinc}(\omega r)$
belongs to a different spectrum and was not the covariance of this sampler.

The implementation filters the FFT of real Gaussian white noise and applies
a deterministic ensemble-variance normalization. It retains the random DC
mode and does not standardize each realisation by its own mean and SD.
Those sample-dependent operations would change the Gaussian ensemble.

On a finite FFT grid, the exact covariance is the inverse discrete spectrum,
with periodicity, truncation, and the implemented removal of even-grid
Nyquist planes. The continuum curve is a reference; agreement is not exact
by construction. Empirical autocorrelation is normalized by the true
zero-lag sample value, not the average of a radial bin.

## Outputs and analytical scope

The corrected script writes `results_corrected/` and refuses a nonempty
destination. Historical figures under `results/` retain the earlier
normalization and reference curve and have not been regenerated.

A valid spectrum alone does not demonstrate higher nonspanning aggregate
VF. A GKF calculation would need the corresponding component spectral
moment and the usual smoothness/regularity assumptions, followed by a
separate connectivity study.

Reference: [Rasmussen & Williams, GPML, covariance functions](https://gaussianprocess.org/gpml/chapters/RW4.pdf).
