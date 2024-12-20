# rheomelt-parameterization

## Automated Assessment of Viscosity, Shear Modulus, and Bulk Modulus Parameters for Kervazo et al. (2021)

This repository contains Python code designed to automate the assessment of parameters required for evaluating the melt-dependence parameterization of viscosity, shear modulus, and bulk modulus as described by Kervazo et al. (2021)
 [https://doi.org/10.1051/0004-6361/202039433]:

$\bullet(\phi) = \bullet_l \frac{1 + \theta^\delta}{(1 - F(\theta, \xi, \gamma))^{B(1 - \phi_\star)}}$

$\theta = \frac{1 - \phi}{1 - \phi_\star}$

$F(\theta, \xi, \gamma) = (1 - \xi) \cdot \text{erf}\left(\frac{\sqrt{\pi}}{2(1 - \xi)}\theta(1 + \theta^\gamma)\right)$

- $\bullet$: Viscosity ($\eta$), shear modulus ($\mu$), or bulk modulus ($K$).
- $\bullet_l$: Viscosity ($\eta_l$), shear modulus ($\mu_l$), or bulk modulus ($K_l$) at the liquidus.
- $\phi$: Melt fraction.
- $\phi_\star$: Parameterization of the transition between solid-like and liquid-like behavior.
- $\delta$: Power-law dependence of the rheological parameters on the melt fraction ($\phi$), describing behavior for melt fractions below the critical melt fraction ($\phi_c$).
- $\gamma$: Transition width between regimes.
- $\xi$: Describes the decrease in strength with increasing melt fraction.
- $B$: Einstein coefficient (set to 2.5).

The parameters $(\phi_\star, \delta, \xi)$ depend on specific rheological parameters, including:

- $\bullet_s$: Viscosity ($\eta_s$), shear modulus ($\mu_s$), or bulk modulus ($K_s$) at the solidus.
- $\phi_c$: Rheological critical melt fraction.
- $a$: Parameter describing the viscosity dependence for $\phi < \phi_c$, i.e., $\eta(\phi) = \eta_s \exp(-a\phi)$.
- $A$: Parameter describing the dependence of shear and bulk moduli on $\phi$ for $\phi < \phi_c$.

    $\bullet(\phi)=\bullet_s(1+A\phi)^{-1}$

    $\mu(\phi)=\mu_s\left(1+\phi\left(\frac{40-24\nu}{15}\right)\right)^{-1}$

    $K(\phi)=K_s\left(1+\phi\left(\frac{5-4\nu}{3(1-2\nu)}\right)\right)^{-1},$

    where $\nu$ is the Poisson ratio (default value $\nu=0.25$).

The code provides an iterative assessment of parameters:

- **Viscosity parameters** ($\delta, \xi, \phi_\star$) derived from ($\eta_s, \phi_c, a$).
- **Shear modulus parameters** ($\delta, \xi, \phi_\star$) derived from ($\mu_s, \phi_c, A$).
- **Bulk modulus parameters** ($\delta, \xi, \phi_\star$) derived from ($K_s, \phi_c, A$).

### Modes

- **one-run mode**: Computes Kervazo-like parameters for a single set of input values specified directly within the code.
  - **usage**: `python3 parameter.py one-run`
- **multi-run mode**: Reads and processes multiple parameter sets from an input file.
  - **Usage**: `python3 parameter.py multi-run <input_file>`
  - **Example** of input file: `data-test.in`
  - **Input File Layout**:
    - The first line serves as a header and is skipped.
    - Starting from the second line, each subsequent line contains parameters that need to be converted
    - Subsequent lines contain:
      - Column 0:  $\eta_l$
      - Column 1:  $\eta_s$
      - Column 2: $a$
      - Column 3: $\mu_l$
      - Column 4: $\mu_s$
      - Column 5: $K_l$
      - Column 6: $K_s$
      - Column 7: $\phi_c$
      - Column 8: $\gamma$
  - **Output Files Layout**
    - `eta-<input_file>`
      - Line 0: header
      - Starting from the second line, each subsequent line contains converted parameters corresponding to the input.
      - Column 0: $\eta_l$ (Pa·s) viscosity of liquidus
      - Column 1: $\delta$ for $\eta$
      - Column 2: $\xi$ for $\eta$
      - Column 3: $\gamma$ for $\eta$
      - Column 4: $\phi_\star$ for $\eta$
    - `mu-<input_file>`
      - Line 0: header
      - Starting from the second line, each subsequent line contains converted parameters corresponding to the input.
      - Column 0: $\mu_l$ (Pa) shear modulus of liquidus
      - Column 1: $\delta$ for $\mu$
      - Column 2: $\xi$ for $\mu$
      - Column 3: $\gamma$ for $\mu$
      - Column 4: $\phi_\star$ for $\mu$
    - `k-<input_file>`
      - Line 0: header
      - Starting from the second line, each subsequent line contains converted parameters corresponding to the input.
      - Column 0: $K_l$ (Pa) bulk modulus of liquidus
      - Column 1: $\delta$ for $K$
      - Column 2: $\xi$ for $K$
      - Column 3: $\gamma$ for $K$
      - Column 4: $\phi_\star$ for $K$


### Computational scheme

- 0-th iteration:
  - i) Initial value of $\phi_\star = {\phi_\star^{\rm ini}}$, with a default value of 0.5.
  - ii) Computation of $\delta$ from $a$ (viscosity) or $A$ (shear and bulk moduli), assuming knowledge of $\phi_\star$ (see below).
  - iii) Evaluation of $\xi$ from the solidus and liquidus values, assuming knowledge of $\phi_\star$ and $\delta$ (see below).
- i-th iteration:
  - i) Estimate of $\phi_\star$ based on $\delta$ and $\xi$ from the previous step (see below).
  - ii) Computation of $\delta$ from $a$ (viscosity) or $A$ (shear and bulk moduli), assuming $\phi_\star$ from the current step (see below).
    iii) Evaluation of $\xi$ from the solidus and liquidus values, assuming knowledge of $\phi_\star$ and $\delta$ from the current step (see below).
- End of iterations:
  - The iterative process stops when the relative change in $\phi_\star$ to be lower than $\epsilon$, i.e.  $\left|\frac{\phi_\star^i-\phi_\star^{i-1}}{\phi_\star^{i-1}}\right|<\epsilon$, default value of $\epsilon$ is set to $10^{-8}$.

### $\phi_\star$ estimate

From the viscosity dependence, we can separate a trend describing the viscosity behavior for low porosities, denoted as $\eta^{LP}$

$\eta^{LP}(\phi)=\eta_l\frac{(1+\theta^\delta)}{\xi^{B(1-\phi_\star)}}.$

We find $\phi_\star$ such that the total trend $\eta(\phi)$ and the partial trend $\eta^{LP}(\phi)$ for low $\phi$ differ by only a fraction of the solidus value at the point $\phi_c$, i.e.,

$\eta(\phi_c)-\eta^{LP}(\phi_c)=f\eta_s,$

where $f$ is a numerical parameter representing the fraction, with a default value of $f=0.1$.

### $\delta$ from $a$ and $\phi_\star$

We know that the viscosity should follow

$\eta(\phi)=\eta_{s}\exp(-a\phi)$

for $\phi<\phi_c$. Therefore, the slope of $\ln\eta(\phi)$ is $-a$. Moreover, for $\phi<\phi_c$, $\textrm{erf}(x)$ is almost constant. We can express the viscosity dependence as follows:

$\eta(\phi)=C(1+\theta^\delta)$.

The easiest way to derive $\delta$  (aside from fitting) is to assume that the two formulas yield the same values for two different values of $\phi$ the two formulas have the same values. The analytical assessment is straightforward for $\phi=0$ and $\phi=\phi_\star$.

The first guess is given by:

$\eta(\phi=0)=\eta_s=C\left(1+\frac{1}{(1-\phi_\star)^\delta}\right),$

$\eta(\phi=\phi_\star)=\eta_s\exp(-a\phi_\star)=2C$

together this leads to

$2\exp(a\phi_\star)=1+\frac{1}{(1-\phi_\star)^\delta}$

$\delta=-\frac{\ln(2\exp(a\phi_\star)-1)}{\ln(1-\phi_\star)}$

As $\phi_\star>\phi_c$, this leads to underestimation of the parameter $\delta$. Therefore, we can apply an iterative improvement for $\phi=\phi_c$ using

$\eta(\phi=0)=\eta_s=C\left(1+\frac{1}{(1-\phi_\star)^\delta}\right)$

$\eta(\phi=\phi_c)=\eta_s\exp(-a\phi_c)=C\left(1+\frac{(1-\phi_c)^\delta}{(1-\phi_\star)^\delta}\right)$

$((1-\phi_\star)^\delta+1)\exp(-a\phi_c)=(1-\phi_\star)^\delta+(1-\phi_c)^\delta.$

Assuming $\delta$ large and $0<1-\phi_\star<1$, we can simplify the equation using $(1-\phi_\star)^\delta+1\approx 1$

$\exp(-a\phi_c)=(1-\phi_\star)^\delta+(1-\phi_c)^\delta.$

$\delta=\frac{\ln(\exp(-a\phi_c)-(1-\phi_\star)^{\delta})}{\ln(1-\phi_c)}$

Iteration scheme is as follow:

- 0th iteration: $\delta_0=-\frac{\ln(2\exp(a\phi_\star)-1)}{\ln(1-\phi_\star)}$
- i-th iteration: $\delta_i=\frac{\ln(\exp(-a\phi_c)-(1-\phi_\star)^{\delta_{i-1}})}{\ln(1-\phi_c)}$.

This iteration scheme converges in just a few steps.

### $\delta$ from $A$ and $\phi_\star$

Let us assume the following interpretation, which is applicable to both moduli:

$f(\phi)=f_s\left(1+A\phi\right)^{-1},$

where $A = \frac{40 - 24\nu}{15}$ for shear modulus and $A = \frac{5 - 4\nu}{3(1 - 2\nu)}$ for bulk modulus, respectively. The notation $f_s$ represents the corresponding solidus values. In a similar fashion to viscosity, we aim to determine $\delta$ for moduli

$g(\phi)=C(1+\theta^\delta)=C\left(1+\left(\frac{1-\phi}{1-\phi_\star}\right)^\delta\right),$

where $g$ denotes the desired representations for low fractions of melt. We can establish the relationship through the minimization problem given by $\min_{{C, \delta}}||f(\phi) - g(\phi)||$ for $\phi \in [0, \phi_c]$.

The most practical approach is to adhere to the conditions $f(\phi=0) = g(\phi=0)$ and $f(\phi=\phi_\star) = g(\phi=\phi_\star)$. Under these conditions, $\delta$ is independent of $f_s$, allowing us to derive an analytic solution.

The condition for $\phi=0$ follows

$f(0)=f_s=C\left(1+\frac{1}{(1-\phi_\star)^\delta}\right)=g(0)$

and $\phi=\phi_\star$ to

$2C=f_s(1+A\phi_\star)^{-1}=C\left(1+\frac{1}{(1-\phi_\star)^\delta}\right)(1+A\phi_\star)^{-1}$

and combined

$1-2\phi_\star A=\frac{1}{(1-\phi_\star)^\delta}\quad\rightarrow\quad
    \delta=-\frac{\ln(1-2A\phi_\star)}{\ln(1-\phi_\star)}.$

However, since $\phi_c < \phi_\star$, the fit may not be satisfactory. We can iteratively derive values of $\delta$ based on the conditions at $\phi = 0$ and $\phi = \phi_c$, specifically:

$f(0)=f_s=C\left(1+\frac{1}{(1-\phi_\star)^\delta}\right)=g(0)$

$f(\phi_c)=f_s(1+A\phi_c)^{-1}
  =C\left(1+\left(\frac{1-\phi_c}{1-\phi_\star}\right)^\delta\right)
    =g(\phi_c)$
$\left(1+\frac{1}{(1-\phi_\star)^\delta}\right)\frac{1}{1+A\phi_c}=\left(1+\frac{(1-\phi_c)^\delta}{(1-\phi_\star)^\delta}\right)$
$\delta=\frac{\ln\left(\frac{1+(1-\phi_\star)^{\delta}}{A\phi_c+1}-(1-\phi_\star)^{\delta}\right)}
    {\ln(1-\phi_c)}$.

The first iteration corresponds to the previous solution. Thus, the iterative scheme is as follows:

- 0th iteration: $\delta_1=-\frac{\ln(1-2A\phi_\star)}{\ln(1-\phi_\star)}$
- i-th iteration: $\delta_i=\frac{\ln\left(\frac{1+(1-\phi_\star)^{\delta_{i-1}}}{A\phi_c+1}-(1-\phi_\star)^{\delta_{i-1}}\right)}
    {\ln(1-\phi_c)}$.

### $\xi$ from $\bullet_s$, $\delta$ and $\phi_\star$

Assuming we already know $\phi_\star$ and $\delta$, we aim to derive $\xi$ from $g_s = g(\phi=0)$, where $g$ represents $\eta$, $\mu$, and $K$. Here, $g_s$ and $g_l$ denote the corresponding solidus and liquidus values, respectively. As $\phi$ approaches 0, we note that $\textrm{erf}(x)$ approaches 1.

$g_s=\eta(\phi=0)=g_l\frac{1+\theta^\delta}{\xi^{B(1-\phi_\star)}}$

$\xi=\left(\frac{g_l}{g_s}(1+\theta(\phi=0))^\delta\right)^{\frac{1}{B(1-\phi_\star)}}.$
