# stochastic-linearalgebra

Research implementation, in MATLAB, of stochastic numerical linear algebra methods. The code in
[`Random Projection/`](Random%20Projection/) builds random projection matrices of Johnson–
Lindenstrauss type and uses them to estimate the `k` entries of largest magnitude of a matrix
product `A·B` from projected inner products of the rows of `A` and the columns of `B`.

## Running

The repository does not record which MATLAB release the code was written for; it dates from
May 2014 and declares no toolboxes. From MATLAB, change into the `Random Projection`
directory and run the test script:

```matlab
cd 'Random Projection'
test_script
```

[`test_script.m`](Random%20Projection/test_script.m) draws random matrices `A` of size
`2000 × 1000` and `B` of size `1000 × 2000` with independent normal entries of standard
deviation `1000`, times `getKLargestMatrixMultiply` against the brute-force computation of
`A*B`, compares their outputs, and times the construction of the projection matrices.

To use the method directly:

```matlab
kLargestElems = getKLargestMatrixMultiply(A, B, k, projection_type)
```

where `projection_type` is one of `'FJLT'` (the default), `'JLT'`, `'ACH'`, `'LHC'` and
`'None'` (a zero matrix, used to measure the baseline error rate).

## Repository layout

| File in `Random Projection/` | Contents |
|---|---|
| [`getKLargestMatrixMultiply.m`](Random%20Projection/getKLargestMatrixMultiply.m) | the `k` largest entries in magnitude of `A·B`, through a chosen projection |
| [`JLT.m`](Random%20Projection/JLT.m) | standard Johnson–Lindenstrauss transform |
| [`ACH.m`](Random%20Projection/ACH.m) | sparse transform of Achlioptas (2001) |
| [`LHC.m`](Random%20Projection/LHC.m) | very sparse transform of Li, Hastie and Church (2006) |
| [`FJLT.m`](Random%20Projection/FJLT.m) | fast Johnson–Lindenstrauss transform of Ailon and Chazelle |
| [`project.m`](Random%20Projection/project.m) | projected inner products of the rows of `A` and the columns of `B` |
| [`test_script.m`](Random%20Projection/test_script.m) | speed and accuracy test |

## Authors and citation

The code is by **Ahmed Bou-Rabee**. There is no associated paper, so the repository has no
`CITATION.cff`.

## License

The code is licensed under the Apache License 2.0; see [`LICENSE`](LICENSE).
