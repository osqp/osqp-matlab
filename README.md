# Matlab interface for OSQP

[![Matlab Interface Tests](https://github.com/oxfordcontrol/osqp-matlab/actions/workflows/ci.yml/badge.svg)](https://github.com/oxfordcontrol/osqp-matlab/actions/workflows/ci.yml)

Matlab wrapper for [OSQP](https://osqp.org/): the Operator Splitting QP Solver.

The OSQP (Operator Splitting Quadratic Program) solver is a numerical optimization package for solving problems in the form
```
minimize        0.5 x' P x + q' x

subject to      l <= A x <= u
```

where `x in R^n` is the optimization variable. The objective function is defined by a positive semidefinite matrix `P in S^n_+` and vector `q in R^n`. The linear constraints are defined by matrix `A in R^{m x n}` and vectors `l in R^m U {-inf}^m`, `u in R^m U {+inf}^m`.


## Documentation
  The installation procedure is documented in [https://osqp.org/docs/get_started/matlab.html](https://osqp.org/docs/get_started/matlab.html), while the interface is documented in [https://osqp.org/docs/interfaces/matlab.html](https://osqp.org/docs/interfaces/matlab.html).
