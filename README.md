# AubryAndre J and Delta from Lattice Parameters

There are four relevant functions:
- __J_Delta_Gaussian__ — the roughest approximation of J and Delta that basically assumes each well is a harmonic oscillator (with gaussian wannier states)
- __J_Numeric_Approx__  —  this is a specific approximate formula for J referenced in [the Modugno paper](https://iopscience.iop.org/article/10.1088/1367-2630/11/3/033023).  It’s better than the previous function by a bit.  The paper it is from claims that it’s accurate to ~1% for primary lattice depths between 8Er and 30Er, but I find that it seems to deviate by 2-3% from the exact numerical calculations, and I’m inclined to trust the numerics, since they converge quite clearly and it was written by someone who would have stress tested I believe.
- __J_Delta_Numeric__ — calculates J and Delta by numerical determining the maximally localized wannier functions using the MLGWS framework (based on Wannier90).  I’ll note that the calculation of J is done by MLGWS in the Bloch basis — it was built into MLGWS.
- __compare_J_Delta_Methods__ — Generates a table that compares the values of the different methods if you put in an array of parameters (s1, lambda1, or lambda2 — s2 is an option, but it only scales Delta)
