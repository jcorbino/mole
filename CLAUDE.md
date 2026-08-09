# MOLE — Mimetic Operators Library Enhanced

High-order **mimetic finite-difference** operators for solving PDEs, in C++ and
MATLAB/Octave from a single shared design. The operators are sparse matrices —
Divergence, Gradient, Laplacian, Bilaplacian, Curl, and interpolators — acting on
**staggered grids** (uniform, non-uniform, curvilinear) and satisfying local and
global conservation laws by construction. Orders of accuracy k = 2, 4, 6, 8.

Mathematics: Corbino & Castillo (2020), <https://doi.org/10.1016/j.cam.2019.06.042>;
earlier matrix-analysis derivation in Castillo & Grone (2003). GPL.

## Layout

| Path | What |
|---|---|
| `mole_MATLAB/` | The MATLAB/Octave library — one function per file, flat |
| `mole_C++/` | The C++ library — operators are classes deriving from `arma::sp_mat` |
| `examples_MATLAB/`, `examples_C++/` | Runnable examples; the C++ ones plot via gnuplot |
| `tests_MATLAB/`, `tests_C++/` | Nullity, energy and accuracy tests |
| `doc_MATLAB/`, `docs/`, `Doxyfile` | Documentation |

The two languages are ports of each other. **A change to one usually needs the
same change to the other** — and the C++ signature often differs, see below.

## The staggered grid, and the one convention to internalise

For an m×n grid, a **scalar** carries `(m+2)(n+2)` values: the m interior cell
centres per axis, plus one entry at each boundary node. A **vector** field is
split by component — `u` on vertical faces is `(m+1) × n`, `v` on horizontal
faces is `m × (n+1)` — so `div2D` maps `(m+1)n + m(n+1) → (m+2)(n+2)`.

**Fields are flattened with x running fastest.** In MATLAB that is exactly what
`ndgrid` produces, so use `ndgrid` and not `meshgrid`, and transpose only when
handing data to a plotting routine. Coordinates:

```matlab
xc = [a a+dx/2 : dx : b-dx/2 b];   % cell space, m+2 long
xf = a : dx : b;                   % faces/nodes, m+1 long
```

## Boundary conditions

Two mutually exclusive styles:

- **Default** — one-sided boundary closures. The boundary rows of `div2D`/`div3D`
  are deliberately *incomplete* (they carry only the normal term, and the corners
  are empty); you are expected to add a BC operator such as `robinBC2D` or
  `mixedBC2D`, which overwrites them.
- **Periodic** — pass `'periodic'` as the trailing `bc` argument to `div`, `grad`,
  `interpol`, their 2-D/3-D counterparts, and `lap`/`lap2D`/`lap3D`. In C++ the
  same flag is a trailing `bool periodic = false`. The interior staggered stencil
  wraps around the seam, no BC operator follows, and every row is a genuine
  operator. MATLAB additionally accepts a **per-axis** form — `{'periodic','none'}`
  in 2-D, three entries in 3-D — for channel and duct geometries; C++ takes a
  single bool for all axes.

With periodic operators the boundary entries of a scalar are duplicates of the
seam: they are carried along for plotting and are never fed back into the
interior. Discrete mass is conserved to roundoff.

## Sharp edges worth knowing before you debug something

- **The interpolators are 2nd-order only at `c = 0.5`.** Any off-centre
  coefficient is 1st-order; upwinding (`c = 1`) and downwinding (`c = 0`) are
  1st-order by design, and that truncation term *is* their numerical diffusion.
- **`interpol` has no order argument at all.** So a composite like
  `div * interpol` is 2nd-order however large `k` is. Raising `k` does not
  improve an advection scheme built that way — it only helps the diffusive term.
- **The default (non-periodic) Laplacian is order k−1**, not k, in the max norm.
  The loss is entirely localised at the rows touching the one-sided closures; the
  deep interior is order k. The periodic Laplacian is uniformly order k.
- **C++ `Divergence` asserts `k < 7`** while `Gradient` goes to k = 8. MATLAB has
  no such gap.

## Building and testing

```
make                      # recurses into mole_C++, examples_C++, tests_C++
make -C tests_C++ run     # C++ tests
```

MATLAB/Octave: `cd tests_MATLAB && run_tests` — Octave runs the library fine, so
`octave-cli` is a valid way to check MATLAB-side changes.

C++ gotchas:

- The Makefile puts `-fopenmp` in `CXXFLAGS`, which **Apple clang rejects**
  without `libomp`. To test on macOS, compile the sources directly rather than
  editing the Makefile.
- `Divergence`/`Gradient`/`Interpol`/`Laplacian` derive from `sp_mat`, but
  Armadillo's expression templates **reject the derived types** — bind to the base
  (`sp_mat D = Dm;`) before multiplying.
- C++ constructors **group** their arguments where MATLAB interleaves them:
  `Divergence(k, m, n, dx, dy)` vs `div2D(k, m, dx, n, dy)`. Easy to transpose by
  accident; the resulting operator is silently wrong rather than an error.

## Conventions

- Commit subjects are imperative and sentence-case, no trailing period, e.g.
  *"Add periodic BC support to the 2-D/3-D operators"*. Most commits are
  subject-only. MATLAB and C++ changes are usually separate commits.
- New arguments go **last and defaulted**, so existing call sites keep working.
- When touching an operator, verify the untouched path is **bit-for-bit**
  unchanged: extract the old file with `git show HEAD:path`, build both, and diff
  the assembled matrices. For C++, dumping sparse triplets at full precision and
  comparing against the MATLAB operator in Octave is what proves a 0-indexed port
  is right.
- Don't commit generated media. Images and PDFs already in `images/` and `docs/`
  are ~12 MB of the repo's ~17 MB.
