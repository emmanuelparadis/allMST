# Find All MSTs Compatible With a Data Set

- [Background](#background)
- [Implementation and Installation](#implementation-and-installation)

## Background

This repository hosts R and C code to find all minimum spanning trees (MST) compatible with a set of distances. ‘Compatible’ means here that all MSTs have the same shortest possible length. The fact that there are several such MSTs is a (potential but not obligate) consequence of ties among the pairwise distances. If all distances are unique, the MST is unique as well.

The function `allMST()` has main argument `d`, a distance matrix, and the optional argument `log = FALSE` which is explained below:

```r
R> source("allMST.R")
R> args(allMST)
function (d, log = FALSE, quiet = FALSE)
NULL	
```

The output is a network as coded by [pegas](https://github.com/emmanuelparadis/pegas/) with an additional attribute `"nMST"` which is the number of MSTs compatible with `d`. As a simple example, let's take the woodmouse data from [ape](https://github.com/emmanuelparadis/ape/):

```r
R> library(ape)
R> data(woodmouse)
R> d <- dist.dna(woodmouse, "n")
R> x <- allMST(d, quiet = TRUE)
R> attr(x, "nMST")
[1] 144
```

We recall that the number of possible spanning trees (i.e., whatever their lengths) for a set of pairwise distances among $n$ observation is $n^{n - 2}$. There are 15 sequences in the woodmouse data:

```r
R> 15^13
[1] 1.946195e+15
```

We now build a data set of 15 points where all distances are equal to one (except on the diagonal of course). We expect all possible MSTs to be compatible with the data:

```r
R> d <- matrix(1, n, n)
R> diag(d) <- 0
R> x <- allMST(d, quiet = TRUE)
R> attr(x, "nMST")
[1] 1.946195e+15
```

With more observations, say $n=200$, even calculating the number of possible spanning trees is not possible with standard computations:

```r
R> n <- 200
R> n^(n - 2)
[1] Inf
```

Instead, we could use, among many possible alternative approaches, the function `LargeNumber()` in [ape](https://github.com/emmanuelparadis/ape/):

```r
R> LargeNumber(n, n - 2)
approximately 4.017345 * 10^455
```

The option `log = TRUE` of `allMST()` counts the number of MSTs on the logarithmic scale:

```r
R> d <- matrix(1, n, n)
R> diag(d) <- 0
R> x <- allMST(d, log = TRUE, quiet = TRUE)
R> attr(x, "nMST")
[1] 1049.067
attr(,"logarithm")
[1] TRUE
R> LargeNumber(exp(1), attr(x, "nMST"))
approximately 4.017345 * 10^455
```

For details on the algorithm, refer to the code. This is work in progress (see some comments in the code).

For another approach on the same problem:

Takeo Yamada, Seiji Kataoka & Kohtaro Watanabe (2010) Listing all the minimum spanning trees in an undirected graph. *International Journal of Computer Mathematics* **87**: 3175-3185. Doi: [10.1080/00207160903329699](https://doi.org/10.1080/00207160903329699).

A manuscript has been submitted with details and more references.

## Implementation and Installation

An implementation is available with the function `rmst()` in the package [pegas](https://github.com/emmanuelparadis/pegas/).

New code are provided in this repository with three versions:

- The R function `allMST()` (in the file 'allMST.R') which works as illustrated in the previous section;
- The R function `bigallMST()` (in the same file 'allMST.R') for big data sets when the full distances cannot fit into computer memory (warning: the code is not optimised and can be very slow but it should work for data sets with around 100,000 nodes);
- The R function `allMST2()` in the file 'allMST_C.R' which is similar to the first one but using C code.

To use the third one, the C files must be compiled with the following commands:

```
R CMD SHLIB dot.c
R CMD SHLIB -c all_mst.c
R CMD SHLIB all_mst.o det.o -o all_mst.so
```

These new code are provided for testing purposes and are under development.
