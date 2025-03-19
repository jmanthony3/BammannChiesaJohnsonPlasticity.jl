# Bammann1990Modeling

```@meta
CurrentModule = BammannChiesaJohnsonPlasticity
```

The following equations are those employed in the [Bammann (1990)](@cite bammannModelingTemperatureStrain1990) paper that are implemented in the `Bammann1990Modeling` type, constructor, and that method of the kernel function, `update` associate with this type.
As such, the internal equations use the same nomenclature for plastic strain rate, $\mathbf{D}^{p} \equiv \dot{\epsilon}^{(p)}$; the second-rank, deviatoric tensors for Cauchy stress, $\mathbf{\sigma}'$ and kinematic hardening, $\mathbf{\alpha}'$ which is an ISV; the scalar isotropic hardening, $\kappa$, and the constants associated with the dynamic and static recovery temperature equations.
Of major importance is that, although not explicitly listed in the publication, the equations for $h$ and $H$ are included in this implementation (c. f. [DYNA3D User Manual (1993)](@cite whirley1993dyna3d)).

```math
\begin{aligned}
    % plastic strain rate
    \mathbf{D}^{p} &= f(\theta)\sinh\left[ \frac{ |\mathbf{\xi}| - \kappa - Y(\theta) }{ V(\theta) } \right]\frac{\mathbf{\xi}'}{|\mathbf{\xi}'|}\text{, let }\mathbf{\xi}' = \mathbf{\sigma}' - \mathbf{\alpha}' \\
    % kinematic hardening
    \dot{\mathbf{\alpha}} &= h\mu(\theta)\mathbf{D}^{p} - [r_{d}(\theta)|\mathbf{D}^{p}| + r_{s}(\theta)]|\mathbf{\alpha}|\mathbf{\alpha} \\
    % isotropic hardening
    \dot{\kappa} &= H\mu(\theta)\mathbf{D}^{p} - [R_{d}(\theta)|\mathbf{D}^{p}| + R_{s}(\theta)]\kappa^{2} \\
    % flow rule
    F &= |\sigma - \alpha| - \kappa - \beta(|\mathbf{D}^{p}|, \theta) \\
    % initial yield stress beta
    \beta(\mathbf{D}^{p}, \theta) &= Y(\theta) + V(\theta)\sinh^{-1}\left(\frac{|\mathbf{D}^{p}|}{f(\theta)}\right) \\
    V(\theta)       &= C_{ 1} \exp\left( -\frac{ C_{ 2} }{ \theta } \right) \\
    Y(\theta)       &= C_{ 3} \exp\left(  \frac{ C_{ 4} }{ \theta } \right) \\
    f(\theta)       &= C_{ 5} \exp\left( -\frac{ C_{ 6} }{ \theta } \right) \\
    r_{d}(\theta)   &= C_{ 7} \exp\left( -\frac{ C_{ 8} }{ \theta } \right) \\
    r_{s}(\theta)   &= C_{ 9} \exp\left( -\frac{ C_{10} }{ \theta } \right) \\
    R_{d}(\theta)   &= C_{11} \exp\left( -\frac{ C_{12} }{ \theta } \right) \\
    R_{s}(\theta)   &= C_{13} \exp\left( -\frac{ C_{14} }{ \theta } \right) \\
    h(\theta)       &= C_{15} \exp\left(  \frac{ C_{16} }{ \theta } \right) \\
    H(\theta)       &= C_{17} \exp\left(  \frac{ C_{18} }{ \theta } \right) \\
\end{aligned}
```

## Types
```@autodocs
Modules = [BammannChiesaJohnsonPlasticity, ContinuumMechanicsBase]
Order   = [:type]
Pages   = ["Bammann1990Modeling.jl"]
```

## Functions
```@autodocs
Modules = [BammannChiesaJohnsonPlasticity, ContinuumMechanicsBase]
Order   = [:function]
Pages   = ["Bammann1990Modeling.jl"]
```

## References
```@bibliography
Pages = []
Canonical = false

bammannModelingTemperatureStrain1990
whirley1993dyna3d
```

## Index
```@index
Modules = [BammannChiesaJohnsonPlasticity, ContinuumMechanicsBase]
Order   = [:type, :function]
Pages   = [@__FILE__]
```
