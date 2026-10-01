# Cross-checks with LoopTools 2.16

Two small Fortran drivers, built against LoopTools 2.16 (www.feynarts.de/looptools):

    gfortran -ffixed-line-length-none -I$LT/include c0_labellings.F -L$LT/lib64 -looptools -o c0
    gfortran -ffixed-line-length-none -I$LT/include d0_points.F    -L$LT/lib64 -looptools -o d0

`c0_labellings.F` evaluates C0(-4,-3,-1; 0,0,1) under all six labellings of its lines, with
both LoopTools versions (a = FF, b = Denner) printed by LoopTools' own compare mode.
Result (2026-10-01, hercules):

| labelling (LoopTools order) | version a (FF) | version b (Denner) |
|---|---|---|
| (-4,-3,-1; 0,0,1) | -1.843 | 2.142 + 11.396i |
| (-1,-3,-4; 0,1,0) | 0.671253105748004 | -0.111 + 2.849i |
| (-3,-1,-4; 1,0,0) | -0.721045180425075 | 0.282 + 3.058i |
| (-4,-1,-3; 0,0,1) | -1.843 | -0.111 + 2.849i |
| (-1,-4,-3; 1,0,0) | 0.671 | 2.142 + 11.396i |
| (-3,-4,-1; 0,1,0) | -0.721 | 2.093 - 4.684i |

The true value is -0.585976809672364722 (feynsage closed form, Mathematica NIntegrate, FeynCalc's
Feynman parametrisation integrated, mpmath). Every invariant is spacelike, so C0 must be real and
negative. Package-X 2.1.1 gives +0.671 and -0.721 for two of the labellings, the same digits as
LoopTools' version a. At nearby points (m2^2 = 0.9801, 1.0201; s2 = -0.9, -1.1) version a agrees
with NIntegrate to 16 digits and version b is still wrong.

`d0_points.F` evaluates the two boxes where Package-X is accurate only to 1e-9 and 1e-8:

| point | feynsage | LoopTools a (FF, default) | LoopTools b (Denner) |
|---|---|---|---|
| (1/4,1,-1/2,0; -1/4,1/4; 2,3/4,11/4,2) | 0.0147906701148701 | 4.857 (wrong) | 0.0147906701148695 |
| (-1/2,-1/4,3/4,0; -1,1/4; 11/4,3,11/4,3/4) | 0.00607769685971634 | 0.0060776968597161 | |
