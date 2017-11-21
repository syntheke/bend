# bend
Bend non-positive matrix to make it positive definite


#### Compile using MKL:

```sh
ifort -c -O3 -g -traceback -i8 -heap-arrays  -openmp bend.f90
```

```sh
ifort  -o bend_b -static -O3 -g -traceback -openmp bend_b.o ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a \
${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,\
--start-group  ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a \
${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_ilp64.a -Wl,--end-group -liomp5\
-lpthread -lm -lz
```

#### Run:
```sh
bend -p fam_file -g grm_file -threads numer_of_threads
```

#### Tell me more:

Optimised version of routine in [MTG2](https://sites.google.com/site/honglee0707/mtg2).

Lee, SH and van der Werf, JHJ (2016) MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Bioinformatics 32, 1420-1422
