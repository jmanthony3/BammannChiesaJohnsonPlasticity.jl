# Primer
Bammann-Chiesa-Johnson (BCJ) plasticity models the inelastic deformation of solid materials at a given temperature and strain rate with Internal State Variables (ISVs).
The BCJ-plasticity model satisfies the _Clausius-Duhem inequality_ by placing thermodynamic constraints on ISVs and their rate equations ([Coleman and Gurtin (1967)](@cite colemanThermodynamicsInternalState1967)) which are motivated by equations for dislocation mechanics.
The model assumes that for any continuum mechanics configuration each material point at each moment in time satisfies some equilibrium criterion for a triplet point composed of the deformation gradient, temperature, and current set of ISVs.
The ISVs in the current configuration intrinsically represent their evolution because of the material history afforded by the rate equations.
These ISVs and their rate equations are used to increment the current configuration from the previous configuration by a trial guess that the material deforms elastically and correcting this trial guess for any plastic deformation accrued in the increment: this is the scheme of the _Radial-Return Method_ from [Krieg and Krieg (1977)](@cite kriegAccuraciesNumericalSolution1977).
Using ISVs to model inelastic deformation is not unique to BCJ-plasticity model.
For a broader review of modeling inelastic deformation with ISVs, as well as, a history of the BCJ-plasticity model, the reader is directed to this review paper: [Horstemeyer and Bammann (2010)](@cite horstemeyerHistoricalReviewInternal2010).

## References
```@bibliography
Pages = []
Canonical = false

colemanThermodynamicsInternalState1967
kriegAccuraciesNumericalSolution1977
horstemeyerHistoricalReviewInternal2010
```