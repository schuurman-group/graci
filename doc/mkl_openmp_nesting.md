# MKL threading inside PySCF's OpenMP regions

**Status:** fixed 2026-08-03 (containment). Underlying mechanism not fully
understood; expect to revisit.

**Affects:** any GRaCI input containing more than one `$scf` section, on a
build linked against threaded MKL. Not specific to `$pbdd`.

---

## Summary

PySCF calls BLAS from inside `#pragma omp parallel` regions throughout its
C layer. MKL is supposed to notice this (`omp_in_parallel()`) and run
serially. **Once a bitci calculation has run in the same process, that
detection stops working**: MKL spawns threads inside an already-parallel
region and the result races.

The consequence is silently corrupted MO integrals for every SCF after the
first, which propagates into the CI as a 5–7 Hartree error.

The fix pins MKL to a single thread for the duration of PySCF's threaded
work and restores it afterwards, so bitci keeps threaded BLAS.

---

## Symptoms

What this looked like before it was understood:

- A `$pbdd` generate job failing its reference-energy check by 7 Hartree,
  on one machine but not another.
- The same input giving different answers on repeat runs at 8 threads,
  and correct answers at 1–2 threads.
- A second `$scf` section failing to converge, where the same calculation
  placed first converges normally.
- Nothing reproducible on small molecules.

The CI symptom is worth recognising, because it looks nothing like an
integral problem. Corrupt integrals shift the first MRCI iteration by
roughly 1e-3 Hartree; `refsel = dynamic` then selects a different
reference space (455 configurations rather than 451); that space no longer
spans the target states, `‖Ψ_R‖` collapses from 0.94 to ~0.5, and the
calculation runs away. The MO truncation and the loss of all two-hole
configurations that follow are consequences, not causes.

---

## Root cause

The pattern, e.g. `pyscf/lib/ao2mo/nr_ao2mo.c`:

```c
#pragma omp parallel default(none) shared(ftrans, fmmm, vout, vin, ...)
{
    double *buf = malloc(...);
#pragma omp for schedule(dynamic)
    for (i = 0; i < nij; i++) {
        (*ftrans)(fmmm, i, vout, vin, buf, &envs);   /* -> dgemm_ */
    }
    free(buf);
}
```

`ftrans` dispatches to routines calling `dgemm_` — a dozen call sites in
that file. The same shape appears in the DFT path, which matters more
because it runs in every SCF iteration:

| file | `omp parallel` regions | BLAS calls |
|---|---|---|
| `dft/nr_numint.c` | 11 | 4 |
| `dft/numint_uniform_grid.c` | 2 | 26 |
| `dft/grid_integrate.c` | 2 | 8 |
| `vhf/optimizer.c` | 1 | 8 |
| `ao2mo/nr_ao2mo.c` | 2 | 12 |

The first SCF and its transform are always correct and bit-reproducible.
Every one after a CI calculation is not.

### Evidence

Measured on stilbene, def2-TZVPD, 8 threads (`c1c1.inp`: two identical C1
calculations on separate `$scf` objects).

Correct value for `sum|eri_mo|` is `3.1035e+04`. Observed for the second
SCF across runs: `6.16e4`, `7.87e4`, `8.18e4`, `8.83e4`, `9.10e4`,
`9.16e4`, `9.26e4` — **a different wrong answer every time**.

Two back-to-back calls with byte-identical inputs in the same process:

```
singlet_a   pass1 3.10353491181511345e+04   pass2 3.10353491181511345e+04   IDENTICAL
singlet_b   pass1 6.91132300518105185e+04   pass2 8.09105897828947054e+04   DIFFERENT
```

The corruption was bracketed by instrumentation and localised precisely:

| stage | verdict |
|---|---|
| inputs (`orbs`, `pymol`) | sound — agree to 1e-10, correct molecule |
| **`df.outcore.general`** | **writes a corrupt tensor** |
| readback (blocked and bulk) | faithful |
| `write_integrals` | faithful |
| Fortran read | faithful |

---

## The fix

`graci/core/libs.py`:

```python
def mkl_set_num_threads(n):
    """Set MKL's thread count, returning the previous value."""
    ...

@contextlib.contextmanager
def mkl_single_thread():
    prev = mkl_set_num_threads(1)
    try:
        yield
    finally:
        if prev is not None:
            mkl_set_num_threads(prev)
```

MKL is reached through `MKL_Set_Num_Threads_Local` looked up in the
process image rather than by loading a library by name — MKL is already
resident (bitci links it, and so do PySCF's C extensions), and the file
name differs between an oneAPI layout and a conda one. Returns `None` if
MKL is absent, in which case there is nothing to pin.

Applied at the two places PySCF does heavy threaded work:

- `graci/core/scf.py` — around `mf.kernel()`
- `graci/core/ao2mo.py` — around `df.outcore.general()`

### Why in code rather than the environment

`MKL_NUM_THREADS=1` in a submit script fixes it equally well, and was
verified to do so. It was rejected because it protects only the runs whose
environment happens to be set correctly, and silently stops protecting
anyone who copies a job script without it. The same reasoning led
`mkl_compat.f90` to pin `MKL_NUM_THREADS=1` inside the overlap code's OMP
region (commit `80f0239`).

The guard restores the previous value on exit, so bitci — which calls MKL
from serial context and benefits from threading — is unaffected.

---

## Performance

Stilbene, def2-TZVPD, 626 basis functions, 8 threads:

| | MKL threaded | MKL pinned | ratio |
|---|---|---|---|
| SCF (3 cycles) | 14.4 s | 19.8 s | **1.38×** |
| AO→MO transform | 5.3 s | 5.2 s | 0.98× |

**The transform is free.** All its BLAS sits inside OpenMP regions where
MKL should have been serial anyway; removing the oversubscription cancels
whatever is lost.

**The SCF costs 38%, and all of it is collateral.** Only the XC numerical
integration (`nr_numint.c`) needs the guard. Pinning for the whole
`mf.kernel()` also de-threads the DF-JK build, which calls `dgemm` from
serial context and was never at risk. There is no clean seam from Python
to separate the two.

In a whole job — `Scf.run_pyscf` is ~45% of wall time on a `$pbdd` point
(168.9 s of 372.7 s) — that is **roughly 17% slower overall**.

If this cost ever becomes a problem, narrowing the guard inside
`mf.kernel()` is where the 38% lives.

---

## What is not understood

**Why MKL's nesting detection fails.** The obvious candidate was
`MKL_DYNAMIC`, which governs whether MKL may reduce its thread count —
including to one inside a parallel region. It is not the answer:

```
[MKLSTATE] before transform singlet_a   MKL_DYNAMIC=1  max_threads=8
[MKLSTATE] before transform singlet_b   MKL_DYNAMIC=1  max_threads=8
```

Nothing is flipping the flag and nothing is leaking a thread-count change.
The trigger (a CI calculation having run) and the fix are both confirmed
by experiment; the mechanism inside MKL is not.

**This is why the fix is containment rather than a cure**, and why it is
likely to be revisited. If it is picked up again, that is where to start:
what does a completed bitci calculation change about the process such that
MKL stops recognising it is inside a parallel region?

**Not reportable upstream as-is.** PySCF alone is deterministic — see
`dftest.py` below — so it takes both parties, and is not a PySCF bug they
would action.

---

## Investigation record

Twelve hypotheses were eliminated before the cause was found. Each was
killed by a test rather than an argument. Recorded so they are not
re-investigated.

| hypothesis | how it died |
|---|---|
| A race in bitci's own OpenMP | `twice.inp`: two C1 calculations from one MO array at 8 threads — all 12 scratch files bit-identical, correct answer |
| `hij_disk.f90` write ordering | That file does not run here. The log's `save_hij_double` is `hbuild_double.f90` (serial); DFT/MRCI(2) uses GVVPT2, which is direct |
| `nbuffer` asymmetry | `nextra` is per irrep, so a flattened C1 gets 33 roots where a C2 got 43. Real, but the laptop has the same 33 roots and converges |
| SCF convergence (`conv_tol_grad`) | GRaCI never sets it, so PySCF derives 1e-4 and orbitals are loose (5e-4 in coefficients, energy to 1e-12). Setting 1e-8 made the CI discrepancy *worse*, 2e-6 → 8.5e-4 |
| Single-precision aliased write in `write_integrals` | Tested in isolation on numpy 2.0.1 and 2.3.5: exact for every shape. NumPy buffers the overlap |
| `precision = single` | `double` did not help |
| MO truncation / `refsel` | `truncate = False` did not help; the reference space was already wrong before any MO drop |
| Leaked `tmp_eri` h5py handle | Fixed anyway as hygiene; corruption survived |
| PySCF's prefetch thread | `lib.misc.ASYNC_IO = False` makes `map_with_prefetch` serial. No effect |
| Dual OpenMP runtimes | The laptop has the mixed GCC+MKL build (both `libgomp` and `libiomp5`) and works; the cluster is pure Intel (`icx`, `libiomp5` only) and fails. `LD_PRELOAD=libiomp5` changed nothing. The correlation runs backwards |
| `max_memory` / block count | GRaCI never passes `max_memory`, so the transform runs at PySCF's `df_outcore_max_memory` default of 2000 MB regardless of the node — `PYSCF_MAX_MEMORY` feeds a different config key and never reaches it. Block count is identical on both machines |
| Disk space in `PYSCF_TMPDIR` | TBs free |
| `MKL_DYNAMIC` | Already 1; never flipped |

### Two things that made it hard

**It looked machine-dependent.** It needs a second `$scf`, enough threads
to lose the race, and a system large enough that the transform threads at
all. Small molecules never show it.

**The damage surfaces far downstream.** A 1e-3 shift in the first MRCI
iteration becomes a 7 Hartree runaway three iterations later, by way of a
discrete change in the selected reference space. Nothing about the failure
points back at the integrals.

### A note on `max_memory`

Not the cause, but worth fixing on its own merits: `df.outcore.general`
takes its `max_memory` default from `__config__.df_outcore_max_memory`
(2000 MB), a *different* config key from `MAX_MEMORY`. `PYSCF_MAX_MEMORY`
never reaches it, so the transformation runs as if it had 2 GB no matter
what the node has. Passing `scf.mol.pymol().max_memory` would let PySCF
size its own blocking sensibly; it applies its own 0.45 derating on top.

---

## Reproducing it

Everything is in `~/projects/graci_testing/pbdd/omp_repro/`.

**`c1c1.inp`** — the reproducer. Two identical C1 calculations on separate
`$scf` objects, def2-TZVPD stilbene, 23 states. With the guard removed and
`OMP_NUM_THREADS=MKL_NUM_THREADS=8`, the second refines to 455
configurations rather than 451 and runs away. With the guard, both give
451 and identical energies.

**`twice.inp`** — the control. Two C1 calculations sharing one `$scf`.
Always correct, because only one transform happens and it happens before
any CI.

**`dftest.py <xyz> [-n N] [--with-bitci]`** — calls `df.outcore.general`
N times on a fixed pseudo-MO matrix. No SCF, no CI; minutes rather than
hours. Returns identical results 3/3 even with `libbitci` loaded, which
is what established that the trigger is a CI having *run*, not bitci being
present.

**`checkfile.py <file>...`** — reads a GRaCI integral file independently
of bitci: dimensions, elements read vs expected, `sum|x|`, `max|x|`,
non-finite count, and first/second-half column sums.

**`fingerprint.sh`, `compare_runs.sh`** — log and bitci-scratch
comparison. Note that bitci rewrites its scratch every DFT/MRCI(2)
iteration, so those files hold only the final state: they show *that* two
runs diverged, not where.

### Diagnostics used and removed

Should this need re-instrumenting, the probes were:

- `[AO2MO-IN]` — orbital and molecule checksums entering the transform
- `[AO2MO]` — `sum|eri_mo|` as PySCF returns it
- `[AO2MO-RD]` — blocked vs bulk readback, to separate a bad file from a
  bad large read
- `[AO2MO-2X]` — the transform run twice with identical inputs; this is
  the one that proved it was a race
- `[READCHK]` (`dep/bitci/src/integrals/df.f90`) — file names,
  dimensions, record structure and checksums as bitci reads them
- `[INTCHK]` (`dep/bitci/src/ci/init/precompute_integrals.f90`) —
  `fii`/`fock`/`Vc`/`Vx` sums and sample raw integrals

`[AO2MO-2X]` is the most informative for the money: it distinguishes a
race from deterministic corruption in a single run.

---

## Follow-ups

- **Sweep for other unguarded PySCF entry points.** The guard covers the
  two heavy ones. Anything else calling into PySCF's C layer after a CI
  has run wants the same treatment — `intor` calls for property
  integrals, the non-DF `ao2mo.incore.full` branch, the interaction and
  spinorbit code. The tell is PySCF C with `#pragma omp parallel` around
  BLAS. Only `vhf`, `dft` and `ao2mo` have been checked.
- **Narrow the SCF guard** if the 17% matters.
- **Pass `max_memory`** to `df.outcore.general` (see above).
- **Re-validate anything computed with multiple `$scf` sections** on a
  threaded-MKL build before this fix.
