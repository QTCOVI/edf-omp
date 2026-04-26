# Quick Start

This page walks through a complete EDF calculation for H₆ — the simplest test case
in the `test/` directory.

## What you need

1. A wavefunction file: `h6.wfn` (GAMESS/Gaussian `.wfn` format)
2. An Atomic Overlap Matrix file: `h6.aom` (from PROMOLDEN or another AOM source)
3. An input file: `h6.edfinp`
4. The `edf-omp` binary in the parent directory

## The input file

```
0
h6.aom
h6.wfn
ngroup 6
1 1
1 2
1 3
1 4
1 5
1 6
end
```

Line by line:

| Line | Content | Meaning |
|------|---------|---------|
| `0` | Record 1 (`ioverlap`) | AOM comes from PROMOLDEN format |
| `h6.aom` | Record 2 (`filedat`) | Name of the AOM file |
| `h6.wfn` | Record 3 (`wfnfile`) | Name of the WFN file |
| `ngroup 6` | Record 4 | Divide the molecule into 6 fragments |
| `1 1` … `1 6` | Records 5.1–5.6 | Each group has 1 basin (atoms 1–6) |
| `end` | Terminator | End of input |

Each `1 N` line means: "group i contains 1 basin, and that basin is atom N."

## Running the calculation

```bash
cd test/
OMP_NUM_THREADS=4 ../edf-omp < h6.edfinp > h6_new.edfout
```

- `OMP_NUM_THREADS=4` uses 4 threads (only effective with `edf-omp`).
- Input is read from **stdin**; output goes to **stdout**.

## What to look for in the output

After the header and wavefunction summary, you will see:

### Spin-split probability table

```
# Spin-Splitted probabilities
# G1(a) G1(b) G2(a) G2(b) ...
#   0.0000000075662346   3  3  0  0  0  0  0  0  0  0  0  0
...
```

Each row: probability followed by alpha and beta electron counts for each group.
`G1(a)` is the number of alpha electrons in group 1, `G1(b)` is beta electrons in group 1.

### Spinless EDF table

After the spin-resolved table, EDF sums over spin to give total electron counts:

```
#   0.0002607553   1  1  1  1  1  1    <- dominant: 1 electron per H
#   0.0001234281   2  1  1  1  1  0    <- one H has 2e, one has 0e
...
#   0.999492  <-- SUM,   56 PROBABILITIES
```

For H₆, the dominant configuration `1 1 1 1 1 1` (one electron per hydrogen) confirms
the delocalized ring character.

### Average populations and delocalization indices

```
<n(1)_alpha>  =  0.499746    <- H1 carries 0.5 alpha electrons
<n(1)_beta>   =  0.499746
...
delta_(2 1)   =  0.44337     <- strong delocalization between adjacent H atoms
delta_(3 1)   =  0.06013     <- weaker across the ring
```

The delocalization index \(\delta(i,j) = 2[\langle n_i n_j \rangle - \langle n_i \rangle \langle n_j \rangle]\) measures electron sharing between fragments.

## Next steps

- [Input Format](input-format.md) — full description of all input records
- [Output Guide](output-guide.md) — detailed explanation of each output section
- [Performance](performance.md) — how to get the most out of `edf-omp` threads
