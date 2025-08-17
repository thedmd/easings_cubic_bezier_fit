# Bezier Curve Fitting Project

Cubic-Bézier approximation of Robert's Penner Easing Functions with more than two digits precision.

This tool does make a legwork and approximate unit Bézier curves control points using:

- Philip J. Schneider's piecewise cubic fitting algorithm ("Graphics Gems" 1990), to make initial guess
- Optionally brute force can be used to piece initial curve control points
- Iterative curve refinement by brute forcing best fit with increasingly smaller window of adjustments

## Results

Program output unit Cubic-Bézier control points in a format `(cp1x, cp1y, cp2x cp2y)`.
> Note Unit Cubic-Bézier expands to `(0, 0), CP1, CP2, (1, 1)`, first and last point is implicit.

Values can be easily plugged to CSS `cubic-bezier(0.404988, -0.000291, 0.746889,  0.494143)`.

Error value is squared difference between approximation and reference function for 100 samples.

```cpp
linear       (0.333333,  0.333333, 0.666667,  0.666667) // error: 0.00000000000000000000
quadIn       (0.404988, -0.000291, 0.746889,  0.494143) // error: 0.00000004035684091910
quadOut      (0.460323,  0.922299, 0.822686,  0.997449) // error: 0.00000154898954060940
quadInOut    (0.476189,  0.035341, 0.523811,  0.964658) // error: 0.00041284620651634067
quadOutIn    (0.476189,  0.917038, 0.523811,  0.082962) // error: 0.00041283086011953049
sinIn        (0.362287, -0.000614, 0.676510,  0.491861) // error: 0.00000004698762095326
sinOut       (0.335369,  0.528073, 0.642641,  0.999652) // error: 0.00000031167956200738
sinInOut     (0.362575, -0.001827, 0.637425,  1.001827) // error: 0.00000026398028387468
sinOutIn     (0.589352,  0.953735, 0.410648,  0.046265) // error: 0.00006070458026140493
expIn        (0.639784,  0.018747, 0.844978, -0.056601) // error: 0.00002687472073739460
expOut       (0.155022,  1.056599, 0.360214,  0.981253) // error: 0.00002687472118436125
expInOut     (0.844114, -0.116792, 0.155886,  1.116792) // error: 0.00567340183663941552
expOutIn     (0.108347,  1.069179, 0.891654, -0.069179) // error: 0.00631927193846550540
circIn       (0.555832,  0.001364, 0.999098,  0.449214) // error: 0.00000055198506271246
circOut      (0.000902,  0.550795, 0.444177,  0.998637) // error: 0.00000055201779274905
circInOut    (0.879116,  0.132567, 0.120884,  0.867433) // error: 0.00887236917283953583
circOutIn    (0.078823,  0.825427, 0.921177,  0.174573) // error: 0.00627826226176189488
cubeIn       (0.333378, -0.000008, 0.666683,  0.000055) // error: 0.00000000002415433220
cubeOut      (0.333350,  1.000055, 0.666712,  0.999991) // error: 0.00000000002431708171
cubeInOut    (0.618693, -0.047936, 0.381307,  1.047936) // error: 0.00021367679521216426
cubeOutIn    (0.333331,  1.000000, 0.666667,  0.000000) // error: 0.00000000000005796166
quarticIn    (0.436496,  0.005838, 0.731336, -0.070799) // error: 0.00000380676780304019
quarticOut   (0.268665,  1.070800, 0.563504,  0.994162) // error: 0.00000380676776716040
quarticInOut (0.708489, -0.096416, 0.291511,  1.096416) // error: 0.00116967797928693371
quarticOutIn (0.243410,  1.048356, 0.756591, -0.048356) // error: 0.00088939334691299533
quinticIn    (0.520802,  0.011918, 0.774258, -0.117217) // error: 0.00002661714888760919
quinticOut   (0.225742,  1.117219, 0.479200,  0.988082) // error: 0.00002661714886089316
quinticInOut (0.769811, -0.128635, 0.230189,  1.128634) // error: 0.00371927255834937198
quinticOutIn (0.182089,  1.080548, 0.817911, -0.080547) // error: 0.00362383145900114188
backIn       (0.333331,  0.000000, 0.666668, -0.567193) // error: 0.00000000000006474762
backOut      (0.333336,  1.567193, 0.666663,  1.000000) // error: 0.00000000000003501250
backInOut    (0.721223, -0.546274, 0.278777,  1.546274) // error: 0.02067101174868633470
backOutIn    (0.250051,  1.310034, 0.749949, -0.310034) // error: 0.00705526661477213525
```

## Bonus

`elasticIn` and `elasticInOut` are split into set of Cubic-Bézier curves.

Start with `LineTo` (`M` in SVG) then series of `CurveTo` (`C` in SVG).
Flip order of points and reverse path can be built.

Each curve has accompanying error of approximation

```cpp
// elastic-in curve segments:
{
    {  0.000000,  0.000000 },
    {  0.000000, -0.000897 }, // error: 0.00000018494911557809
    {  0.007781, -0.000259 },
    {  0.025000,  0.000000 },
    {  0.070783,  0.001105 }, // error: 0.00000005197258090561
    {  0.123203,  0.003963 },
    {  0.175000,  0.000000 },
    {  0.278433, -0.009199 }, // error: 0.00000119017988708947
    {  0.287967, -0.006129 },
    {  0.325000,  0.000000 },
    {  0.425176,  0.024895 }, // error: 0.00000807064407126745
    {  0.435771,  0.018783 },
    {  0.475000,  0.000000 },
    {  0.541920, -0.040144 }, // error: 0.00000327888778883789
    {  0.571872, -0.081376 },
    {  0.625000,  0.000000 },
    {  0.687349,  0.102489 }, // error: 0.00000769037433201447
    {  0.720941,  0.237567 },
    {  0.775000,  0.000000 },
    {  0.832459, -0.256431 }, // error: 0.00000618931062490446
    {  0.870057, -0.693085 },
    {  0.925000,  0.000000 },
    {  0.954601,  0.373762 }, // error: 0.00000944733892538352
    {  0.975766,  0.827874 },
    {  1.000000,  1.000000 },
}

// Fitting elastic-in-out (first half) curve segments:
{
    {  0.000000,  0.000000 },
    {  0.111327,  0.002458 }, // error: 0.00000011547778880529
    {  0.155256,  0.003351 },
    {  0.212500,  0.000000 },
    {  0.353944, -0.011087 }, // error: 0.00000284385117008691
    {  0.371161, -0.016770 },
    {  0.437500,  0.000000 },
    {  0.557508,  0.040354 }, // error: 0.00002027147771412858
    {  0.591539,  0.089201 },
    {  0.662500,  0.000000 },
    {  0.755518, -0.124011 }, // error: 0.00000013867749099593
    {  0.815255, -0.464200 },
    {  0.887500,  0.000000 },
    {  0.931458,  0.280987 }, // error: 0.00000003527866548364
    {  0.966584,  0.768286 },
    {  1.000000,  1.000000 },
}
```

## Building

CMake project is provided for convenance. Alternatively code could be compiled
as one-liner.

### Requirements

- C++20 compatible compiler (GCC 10+, Clang 10+, MSVC 2019+)
- OpenMP support (optional, for parallel optimization)

### Compilation

```bash
# Basic compilation
g++ -std=c++20 -O3 src/main.cpp FitCurves.c -o bezier_fit

# With OpenMP support
g++ -std=c++20 -O3 -fopenmp src/main.cpp FitCurves.c -o bezier_fit

# Visual Studio
cl /std:c++20 /O2 /openmp src/main.cpp FitCurves.c
```

## Mathematical Background

Cubic Bézier curves are defined by four control points P₀, P₁, P₂, P₃:

```text
B(t) = (1-t)³P₀ + 3(1-t)²tP₁ + 3(1-t)t²P₂ + t³P₃
```

For easing curves:

- P₀ = (0, 0) - start point
- P₃ = (1, 1) - end point
- P₁, P₂ - control points that define the curve shape

## References

1. Schneider, Philip J. "An Algorithm for Automatically Fitting Digitized Curves." Graphics Gems, Academic Press, 1990.
2. Penner, Robert. "Easing Equations." [robertpenner.com/easing](http://robertpenner.com/easing/)
3. CSS Cubic Bezier reference implementations
4. [Animation easing functions](https://easings.net/)
5. [CSS Easing Animation Tool](https://matthewlein.com/tools/ceaser)
6. [Cubic Bezier Approximations for Robert Penner Easing Equations](https://github.com/zz85/cubic-bezier-approximations/tree/gh-pages)

