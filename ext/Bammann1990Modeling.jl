export parameters

"""
Constants for temperature equations from [Bammann (1990)](@cite bammannModelingTemperatureStrain1990).
Note: though not explicitly listed in paper, temperature equations `h = C₁₅ * exp(C₁₆ / θ)` and `H = C₁₇ * exp(C₁₈ / θ)` are included (c. f. [DYNA3D User Manual (1993)](@cite whirley1993dyna3d)).
"""
ContinuumMechanicsBase.parameters(::Bammann1990Modeling) = (
    :C₁,    :C₂,    # V
    :C₃,    :C₄,    # Y
    :C₅,    :C₆,    # f
    :C₇,    :C₈,    # r_d
    :C₉,    :C₁₀,   # r_s
    :C₁₁,   :C₁₂,   # R_d
    :C₁₃,   :C₁₄,   # R_s
    :C₁₅,   :C₁₆,   # h
    :C₁₇,   :C₁₈    # H
)