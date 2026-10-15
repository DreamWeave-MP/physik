# physik 0.3.2 benchmarks

The move onto l3i 2.0, measured as an A/B with `bench/bench.jsl`. The numbers in `docs/README.md` and `docs/PROFILE.md` are older, from another machine and another Jess, and are not comparable with these.

## On l3i 2.0

Measured on 2026-10-12, at the end of the run started on 2026-10-11, on an Intel(R) Core(TM) i7-10870H CPU @ 2.20GHz (16 threads), with Jess 0.1.0 on l3i 2.0.2, a verified release build, by

```sh
jess run --package jpk/physik --package-path jpk jpk/physik/bench/bench.jsl '!gets wrong'
```

Base is commit `472f61b`: 0.3.2 as it was, with the suite and the bench added and nothing else changed. Head is `c975c9f`, the last change to the code. Same binary, four rounds inside one hold of the machine's bench lock, base first in the odd rounds and head first in the even ones; the figures are the middle of the four (the mean of the two in the middle), each of which is the bench's own median of seven batches. This is the second such A/B: the first, on `819f278` an hour earlier with more else running, read up to a tenth higher on the rows that moved, and 1.01 to 1.08 on the rows nothing touched. Every workload checks its result against what 0.3.2 computes, to the bit: the contacts every step returned, the candidates, manifolds, warm starts, sleepers and islands of the last one, and a checksum of every body's position, rotation and velocities; for a query, the hits and a checksum of what was hit and where. Base and head pass the same checks, so the table compares two ways of doing the same arithmetic in the same order.

Microseconds per call. A call is what the name says: a world built and stepped that many times, one batch of queries, one drain, one shape built.

| Workload | Base | Head | Base / head |
| --- | ---: | ---: | ---: |
| step: 1,000 spheres piled in the grid, 60 steps | 525,324.6 | 370,996.2 | 1.42 |
| step: 1,000 spheres piled in the grid with contact events, 60 steps | 552,291.2 | 394,141.5 | 1.40 |
| step: 250 boxes tumbling wider than a cell, 30 steps | 236,699.2 | 190,371.4 | 1.24 |
| step: 250 boxes tumbling in the grid, 30 steps | 237,108.2 | 183,837.4 | 1.29 |
| step: 600 spheres, capsules and boxes swept and pruned, 30 steps | 243,119.9 | 164,327.3 | 1.48 |
| step: 2,000 spheres swept and pruned from scattered, 1 step | 23,133.3 | 3,642.5 | 6.35 |
| step: 256 compounds of three pieces, 30 steps | 366,707.6 | 287,730.4 | 1.27 |
| step: 600 capsules and boxes on a heightfield 65 by 65, 30 steps | 222,219.0 | 190,772.8 | 1.16 |
| step: 600 capsules and boxes on a mesh of 8,192 triangles, 30 steps | 384,805.0 | 335,729.8 | 1.15 |
| step: 600 spheres, capsules and boxes on a heightfield 65 by 65, swept and pruned, 30 steps | 215,458.6 | 175,702.7 | 1.23 |
| step: 600 spheres, capsules and boxes on a mesh of 8,192 triangles, swept and pruned, 30 steps | 349,392.0 | 295,507.7 | 1.18 |
| step: 5,000 spheres drifting far apart with contact events, drained, 30 steps | 113,293.9 | 86,923.8 | 1.30 |
| step: 1,000 spheres asleep on the floor, 1 step | 1,634.0 | 1,174.4 | 1.39 |
| step: 60 chains of 8 links on distance, ball and hinge joints, 30 steps | 199,899.8 | 175,957.9 | 1.14 |
| step: 200 swept spheres fired at walls, spheres, compounds and terrain, swept and pruned, 6 steps | 125,664.9 | 120,418.2 | 1.04 |
| step: 200 swept spheres fired at walls, spheres and compounds in the grid, 6 steps | 69,665.8 | 60,804.8 | 1.15 |
| query: 2,048 rays in a batch, as columns | 56,227.1 | 38,679.4 | 1.45 |
| query: 2,048 rays in a batch, as columns, unbounded and downward | 42,851.1 | 31,648.6 | 1.35 |
| query: 2,048 rays in a batch, as records | 54,694.5 | 37,480.0 | 1.46 |
| query: 512 rays one by one | 14,005.9 | 9,501.6 | 1.47 |
| query: 1,024 sphere casts in a batch, as columns | 17,938.5 | 13,300.6 | 1.35 |
| query: 2,048 boxes of a few cells | 26,992.7 | 17,772.3 | 1.52 |
| query: 64 boxes of more cells than the grid is walked for | 37,604.3 | 22,717.2 | 1.66 |
| drain: transforms of 1,001 bodies | 56.8 | 42.2 | 1.35 |
| drain: linear state of 1,000 dynamic bodies, with user data and without handles | 70.4 | 62.6 | 1.12 |
| drain: contacts of the pile, as columns | 88.6 | 89.5 | 0.99 |
| drain: contacts of the pile, as records | 184.5 | 188.9 | 0.98 |
| drain: contact events of the pile, without the ones that stay | 61.5 | 62.5 | 0.98 |
| drain: stats of the pile | 31.5 | 19.2 | 1.64 |
| build: a mesh of 8,192 triangles from tables | 2,726.3 | 2,743.3 | 0.99 |
| build: a mesh of 8,192 triangles from buffers | 2,266.0 | 2,267.9 | 1.00 |
| build: a mesh of 1,408 triangles without indices | 584.5 | 549.5 | 1.06 |
| build: a heightfield 65 by 65 from a table | 54.4 | 39.2 | 1.39 |
| build: a heightfield 65 by 65 from a buffer | 17.6 | 10.1 | 1.73 |
| build: add 5,000 bodies of every primitive, then remove every other one | 8,695.3 | 8,806.2 | 0.99 |
| sync: move the origin under 5,000 bodies, and back | 784.8 | 623.3 | 1.26 |
| sync: apply the transforms of 5,000 bodies | 2,404.9 | 2,313.5 | 1.04 |

**Read the floor first.** The rows nothing touched (contacts drained as columns and as records, contact events, adding and removing bodies, applying transforms) read 0.98 to 1.04 in this run, and a row's four rounds spread up to about a tenth either side of its figure. That is the noise: other agents' tests were running on the machine. A row inside it is "no change": those, the swept spheres against terrain (1.04), and the three mesh builds (0.99, 1.00 and 1.06).

**What moved, and why.** In the order the commits were made; the ratios in passing are from three-round A/Bs taken as each landed.

- *Piles and sweeps, 1.09x to 1.38x at once: the pair filter (`b1e983d`).* A broadphase offers a pair, and most pairs it offers are apart: 21,000 offers a step for the 1,000 spheres, of which 3,500 go on to the narrow phase. 0.3.2 checked liveness, a joint table, body types and masks, and read both shapes, before it compared any bounds. It now compares the bounds first, from a column of sphere radii the world keeps beside its shapes, and reads nothing else for a pair that is apart. The pairs that survive are the same ones in the same order.
- *Piles, another tenth or so (`071ce44`).* Contact preparation and the solver read forty-odd columns from the world for every contact in every iteration. They now bind each column once, with record patterns, and the cache match takes the anchors preparation has just made in place of making them again.
- *The 5,000 drifting spheres, 1.03x to 1.18x (`34c32e9`), then to 1.3x (`ba063aa`).* A body crossing a cell rebuilds the grid's occupied-neighbour pairs: thirteen lookups a cell, nearly all of them misses. A byte for each residue of a packed cell now says a cell is empty without the lookup. The rest is the same binding of columns in the integrator, the island builder, the sleep passes, the contact cache and the drains, which is also where the transform and linear-state drains (1.35x, 1.12x), the origin shift (1.26x) and the sleeping floor (1.39x) come from.
- *Rays, sphere casts and box queries, 1.35x to 1.66x (`52e2bb1`).* A ray walked the grid by probing all twenty-seven cells around each cell it entered, eighteen of which it had probed from the cell before. It now probes the whole neighbourhood once and after that only the cells a step brings in; a body is still tested once, at the same point in the walk, against the same best distance. Box queries pack a row's cell keys by adding one, ask the filter before a lookup, and drop a set of marks that could never fire, a body being in one cell.
- *A sweep that starts disordered, 6.3x (`640ec2e`).* The sweep-and-prune order is insertion-sorted, which is right when motion is coherent and quadratic when 2,000 bodies have just been added. Once the sort has moved more than eight times the bodies there are, what is left is handed to `data.argsort` over a buffer of keys, filled by a buffer sink and applied by a comprehension over a typed slice of the permutation. Both sorts are stable, so they end in the same order; keys that are not numbers, which the two would place differently, keep the insertion sort.
- *Compounds, 1.17x to 1.27x (`ba0c230`).* Two three-piece compounds are nine piece pairs, and each pair placed both its pieces in the world again. A side's pieces are now placed once.
- *Heightfields, 1.4x from a table and 1.7x from a buffer (`c32b927`, `819f278`).* The least and greatest sample are one `min` and one `max` on a typed receiver. That took two tries: see below.
- *`stats`, 1.64x (`99d9002`).* Two `#[...]` counts over columns bound once, where two loops read the world for every body and contact.

**What the dialect and the data plane were good for here, and what they were not.** physik keeps its hot state as Luau vectors in tables, on purpose, and exposes two of those tables to trusted callers, so there is no buffer for the data plane to run over in a step, and none of the step rows above is a data-plane row. They moved because of what is done, not how it is spelled: record patterns made binding a world's columns once a line instead of a paragraph, and that is the dialect's whole share of them. Where there was a buffer, the results are mixed and worth knowing:

- `data.argsort` is the one dramatic row, and it is dramatic only where the insertion sort was wrong for the input.
- The filtered receivers lost. `minLt` and `maxGt` say exactly what the loops said (a sample that is not a number is in neither bound), and over 4,225 samples they took 17.7 µs against the loop's 17.9; over a mesh's packed vertices, six strided passes took 72.7 µs against one loop's 57.4. The unfiltered `min` and `max` take 8.5 µs and 48.6. So the bounds use those, and ask again with the filtered ones only when the answer is not a number, which is what the unfiltered ones give when such a value comes first. Same bounds, bit for bit, pinned by tests that put the value first, last and everywhere.
- A check and a sink lost to one loop. Filling a heightfield's buffer from a table as `assert(all[...])` then `into(buffer, 'f32')[...]` read 0.79 against the loop that checks and writes in one pass, and was taken back out before it was committed.
- A sink where nothing is checked is fine: a mesh without indices numbers its own through `into(triangles, 'u32')[for index in range(...) => index]`, at no measurable change.

Left as loops, deliberately: the drains and the batch queries, which fill several columns in one pass and would need a pass apiece as comprehensions; `compound`, which validates as it builds four lists; the solver, the narrow phase and the sweeps, which are arithmetic with early exits; and `contacts`, already a presized loop, whose records carry the public fields `a` and `b`: single letters the style flags in any line that is written again.

## What 0.3.2 gets wrong

Three workloads run only against head, because 0.3.2 cannot be held to what they compute. 0.3.2 stopped rewriting a sphere's bounds every step under the grid broadphase, and three readers of those bounds were not told:

- the narrow phase against a mesh or a heightfield finds the triangles under a body from its bounds, so a sphere that rolled away from where it was added fell through the terrain (`15e136e`);
- a ray tested a sphere's stale bounds before the sphere, and missed one that had moved (`a5cd514`);
- with more than `directLargeThreshold` oversized bodies, the sweep that pairs them sorted and pruned on stale bounds, and a small sphere passed through twenty walls without a contact (`7df82ea`).

Each is fixed where the bounds are read, with a test that fails on 0.3.2. None of the 37 workloads above changes its result, which is why they were chosen: their terrain rows carry no spheres under the grid, and the rows with spheres on terrain run swept and pruned, where every body's bounds are written every step. These three are the ones 0.3.2 computes differently:

| Workload | Head |
| --- | ---: |
| step: 600 spheres, capsules and boxes on a heightfield 65 by 65, 30 steps, which 0.3.2 gets wrong | 170,828.0 |
| step: 600 spheres, capsules and boxes on a mesh of 8,192 triangles, 30 steps, which 0.3.2 gets wrong | 289,618.3 |
| step: 200 swept spheres fired at walls, spheres, compounds and terrain, 6 steps, which 0.3.2 gets wrong | 120,569.3 |

The swept spheres compute, under the grid, exactly what they compute swept and pruned: the same contacts and the same checksum.

A fourth defect was fixed on the way (`c62e8c7`): a body refused for its rotation had already taken a slot, and left a live body with no rotation behind for the next step to trip over.

## Which modules are native, and why

Jess compiles to machine code only what is marked (`--!native` on a module, `@native` on one
function), and every load pays for what is marked. A module carries `--!native` when a row below
spends its time there and the module is small (12 KB or less); a larger module carries `@native`
on the functions the row spends its time in, and nothing else.

Measured 2026-10-15 with `jess 0.1.0`, the verified release build of jess `07c38cb` on l3i
2.0.3, on the i7-10870H, the one-minute load 1.24 and 1.40 at the start of the two
sessions. In each session the bench ran alternately on this package and on a copy of it with
every module marked, three fresh runs of `bench/bench.jsl` each, every run's row the median of
its batches; a column is the median of the three, and a ratio is against the all-marked runs of
the same session. Which functions a row spends its time in comes from `jess run --profile` of
that row alone, unmarked. The load cost is a fresh `jess run --package <this package>
--package-path jpk` of a script that only requires `@dream/physik`, the two sides alternating, the
median of five each after one warm-up.

| Requiring `@dream/physik` | All marked | This state |
| --- | ---: | ---: |
| Unmarked | 0.117 s | 0.058 s |
| Final | 0.124 s | 0.082 s |

| Row | All marked (ms) | Unmarked / all marked | Final (ms) | Final / all marked |
| --- | ---: | ---: | ---: | ---: |
| step: 1,000 spheres piled in the grid, 60 steps | 342.0 | **2.266×** | 470.8 | **1.377×** |
| step: 1,000 spheres piled in the grid with contact events, 60 steps | 362.3 | **2.286×** | 506.4 | **1.398×** |
| step: 250 boxes tumbling wider than a cell, 30 steps | 175.1 | **2.243×** | 222.5 | **1.271×** |
| step: 250 boxes tumbling in the grid, 30 steps | 167.7 | **2.424×** | 209.1 | **1.247×** |
| step: 600 spheres, capsules and boxes swept and pruned, 30 steps | 149.7 | **1.915×** | 178.8 | **1.194×** |
| step: 2,000 spheres swept and pruned from scattered, 1 step | 4.203 | **1.779×** | 5.289 | **1.258×** |
| step: 256 compounds of three pieces, 30 steps | 265.6 | **2.234×** | 380.2 | **1.431×** |
| step: 600 capsules and boxes on a heightfield 65 by 65, 30 steps | 176.5 | **2.628×** | 264.6 | **1.500×** |
| step: 600 capsules and boxes on a mesh of 8,192 triangles, 30 steps | 326.1 | **2.494×** | 393.1 | **1.205×** |
| step: 600 spheres, capsules and boxes on a heightfield 65 by 65, swept and pruned, 30 steps | 169.7 | **2.342×** | 233.1 | **1.373×** |
| step: 600 spheres, capsules and boxes on a mesh of 8,192 triangles, swept and pruned, 30 steps | 283.0 | **2.452×** | 349.0 | **1.233×** |
| step: 600 spheres, capsules and boxes on a heightfield 65 by 65, 30 steps, which 0.3.2 gets wrong | 178.0 | **2.336×** | 230.6 | **1.295×** |
| step: 600 spheres, capsules and boxes on a mesh of 8,192 triangles, 30 steps, which 0.3.2 gets wrong | 281.4 | **2.422×** | 344.0 | **1.222×** |
| step: 5,000 spheres drifting far apart with contact events, drained, 30 steps | 83.35 | **2.158×** | 86.90 | 1.043× |
| step: 1,000 spheres asleep on the floor, 1 step | 1.301 | **2.228×** | 1.629 | **1.252×** |
| step: 60 chains of 8 links on distance, ball and hinge joints, 30 steps | 168.3 | **2.186×** | 222.9 | **1.324×** |
| step: 200 swept spheres fired at walls, spheres, compounds and terrain, swept and pruned, 6 steps | 113.3 | **2.024×** | 125.1 | **1.104×** |
| step: 200 swept spheres fired at walls, spheres and compounds in the grid, 6 steps | 56.79 | **2.022×** | 62.92 | **1.108×** |
| step: 200 swept spheres fired at walls, spheres, compounds and terrain, 6 steps, which 0.3.2 gets wrong | 114.5 | **2.032×** | 137.3 | **1.199×** |
| query: 2,048 rays in a batch, as columns | 36.45 | **2.618×** | 37.42 | 1.027× |
| query: 2,048 rays in a batch, as columns, unbounded and downward | 29.93 | **2.498×** | 31.49 | **1.052×** |
| query: 2,048 rays in a batch, as records | 35.35 | **2.689×** | 35.96 | 1.017× |
| query: 512 rays one by one | 9.167 | **2.708×** | 9.136 | 0.997× |
| query: 1,024 sphere casts in a batch, as columns | 12.68 | **1.972×** | 13.95 | **1.100×** |
| query: 2,048 boxes of a few cells | 17.90 | **2.053×** | 17.86 | 0.997× |
| query: 64 boxes of more cells than the grid is walked for | 24.76 | **2.220×** | 24.53 | 0.991× |
| drain: transforms of 1,001 bodies | 0.2059 | **1.227×** | 0.1988 | 0.965× |
| drain: linear state of 1,000 dynamic bodies, with user data and without handles | 0.1303 | **1.279×** | 0.1282 | 0.984× |
| drain: contacts of the pile, as columns | 0.2156 | **1.565×** | 0.2246 | 1.042× |
| drain: contacts of the pile, as records | 0.2520 | **1.575×** | 0.2330 | 0.925× |
| drain: contact events of the pile, without the ones that stay | 0.05996 | **3.096×** | 0.05926 | 0.988× |
| drain: stats of the pile | 0.01843 | **2.309×** | 0.01804 | 0.979× |
| build: a mesh of 8,192 triangles from tables | 3.126 | **3.397×** | 3.022 | 0.967× |
| build: a mesh of 8,192 triangles from buffers | 2.759 | **3.059×** | 2.656 | 0.963× |
| build: a mesh of 1,408 triangles without indices | 0.6091 | **3.166×** | 0.5870 | 0.964× |
| build: a heightfield 65 by 65 from a table | 0.03816 | **4.891×** | 0.03725 | 0.976× |
| build: a heightfield 65 by 65 from a buffer | 0.01000 | **1.467×** | 0.00980 | 0.981× |
| build: add 5,000 bodies of every primitive, then remove every other one | 12.46 | **1.495×** | 12.67 | 1.017× |
| sync: move the origin under 5,000 bodies, and back | 1.381 | **1.317×** | 1.349 | 0.977× |
| sync: apply the transforms of 5,000 bodies | 3.073 | **1.964×** | 2.951 | 0.960× |

What carries the directive: `physik`, whole, the package's only module, as before this pass.

Marking only the functions the rows' profiles named (the narrow phase, the solvers, the
integrator, the grid and sweep-and-prune, the casts) saved 0.041 s a load and left the steps
1.10x to 1.50x slower than all marked: 1,000 spheres piled in the grid went from 0.34 s to 0.47 s
for 60 steps. A world steps every frame and loads once, so the module stays marked whole. The
table's final column is that function-marked session, kept as the record of what was tried; the
package ships as all marked.
