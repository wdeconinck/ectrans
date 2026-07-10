# NPROMATR and bit reproducibility

`NPROMATR` controls field blocking in the transform package. It is set through
the public `SETUP_TRANS0` optional argument `KPROMATR`; `NPROMATR=0` means that
all fields are transformed together, while a positive value splits the fields
into packets of at most that size.

Different `NPROMATR` values are mathematically equivalent transform settings,
but they are not guaranteed to be bit-identical. Checksum or bitwise reference
tests must therefore compare results produced with the same `NPROMATR` value.

## Compared benchmark tests

The investigation compared two otherwise identical direct-transform benchmark
runs:

| Setting | Full-field run | `NPROMATR=20` run |
| --- | --- | --- |
| Executable | `ectrans-benchmark-cpu-dp` | `ectrans-benchmark-cpu-dp` |
| Truncation | `--truncation 47` | `--truncation 47` |
| Grid | `--grid O48` | `--grid O48` |
| Call mode | `--callmode 1` | `--callmode 1` |
| Fields and levels | `--nfld 10 --nlev 20` | `--nfld 10 --nlev 20` |
| Iterations | `--niter 1 --niter-warmup 0` | `--niter 1 --niter-warmup 0` |
| Checksum frequency | `--check 100` | `--check 100` |
| Field blocking | default `NPROMATR=0` | `--npromatr 20` |
| Direct transform shape at `KM=0` | full field set, `KIFC=241` | tail packet, `KIFC=1` |

The only intentional difference was the transform field packet size. The traced
nonmatching field was the final scalar field in the `NPROMATR=20` tail packet,
which corresponds to scalar field 201 in the full-field run.

## Why the bits can differ

Changing `NPROMATR` changes the field batching used by the direct and inverse
transform control routines. In the direct transform, `DIR_TRANS_CTL` uses
`FIELD_SPLIT` when `NPROMATR > 0` and then calls the Fourier and Legendre steps
packet by packet. Without splitting, the same transform is called once with the
full field set.

This changes the shapes of the dense linear algebra operations inside the
Legendre transform. A traced `T47 O48` direct transform with `NPROMATR=20` showed
the first checksum difference in `LEDIR` for `KM=0`.

The Fourier input column and the weighted `ZB` input to the Legendre GEMM were
bit-identical. The first difference appeared in the GEMM result `ZCA`. In the
full path the call transformed all active fields together, for example with
`KIFC=241`; in the `NPROMATR=20` tail packet the same scalar field was transformed
alone, with `KIFC=1`.

BLAS implementations are allowed to choose different kernels, blocking, vector
widths, and accumulation orders for different matrix shapes. Floating-point
addition and multiplication are deterministic for a fixed operation sequence,
but changing the operation sequence can change the final rounded bits. The
results should remain within normal floating-point roundoff expectations, but
bit-identical output across different `NPROMATR` values should not be expected.

## Practical guidance

Use the same `KPROMATR` value for any run that must be compared bit-for-bit with
another run. This includes checksum reference generation, regression tests, and
local/external reference comparisons.

Use `KPROMATR` as a performance and memory-layout tuning parameter, not as a
bit-preserving switch. If a reference is needed for `NPROMATR > 0`, generate and
store a reference with that same `NPROMATR` value.