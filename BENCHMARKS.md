# physik 0.3.2 benchmarks

The move onto l3i 2.0, measured as an A/B with `bench/bench.luau`. The numbers in `docs/README.md` and `docs/PROFILE.md` are older, from another machine and another Jess, and are not comparable with these.

## On l3i 2.0

Measured on 2026-10-12, at the end of the run started on 2026-10-11, on an Intel(R) Core(TM) i7-10870H CPU @ 2.20GHz (16 threads), with Jess 0.1.0 on l3i 2.0.1, a verified release build, by

```sh
jess run --package jpk/physik --package-path jpk jpk/physik/bench/bench.luau '!gets wrong'
```

Base is commit `472f61b`: 0.3.2 as it was, with the suite and the bench added and nothing else changed. Head is `819f278`. Same binary, four rounds inside one hold of the machine's bench lock, base first in the odd rounds and head first in the even ones; the figures are the median of the four, each of which is the bench's own median of seven batches. Every workload checks its result against what 0.3.2 computes, to the bit: the contacts every step returned, the candidates, manifolds, warm starts, sleepers and islands of the last one, and a checksum of every body's position, rotation and velocities; for a query, the hits and a checksum of what was hit and where. Base and head pass the same checks, so the table compares two ways of doing the same arithmetic in the same order.

Microseconds per call. A call is what the name says: a world built and stepped that many times, one batch of queries, one drain, one shape built.

| Workload | Base | Head | Base / head |
| --- | ---: | ---: | ---: |
| step: 1,000 spheres piled in the grid, 60 steps | 558,011.8 | 386,126.6 | 1.45 |
| step: 1,000 spheres piled in the grid with contact events, 60 steps | 590,409.5 | 396,842.8 | 1.49 |
| step: 250 boxes tumbling wider than a cell, 30 steps | 252,279.7 | 192,168.1 | 1.31 |
| step: 250 boxes tumbling in the grid, 30 steps | 249,657.9 | 182,937.0 | 1.36 |
| step: 600 spheres, capsules and boxes swept and pruned, 30 steps | 254,141.8 | 163,127.3 | 1.56 |
| step: 2,000 spheres swept and pruned from scattered, 1 step | 24,010.8 | 3,730.0 | 6.44 |
| step: 256 compounds of three pieces, 30 steps | 372,364.9 | 296,567.6 | 1.26 |
| step: 600 capsules and boxes on a heightfield 65 by 65, 30 steps | 230,857.1 | 195,699.1 | 1.18 |
| step: 600 capsules and boxes on a mesh of 8,192 triangles, 30 steps | 399,774.0 | 347,341.1 | 1.15 |
| step: 600 spheres, capsules and boxes on a heightfield 65 by 65, swept and pruned, 30 steps | 222,217.3 | 177,467.3 | 1.25 |
| step: 600 spheres, capsules and boxes on a mesh of 8,192 triangles, swept and pruned, 30 steps | 356,982.0 | 305,231.5 | 1.17 |
| step: 5,000 spheres drifting far apart with contact events, drained, 30 steps | 117,715.4 | 90,540.9 | 1.30 |
| step: 1,000 spheres asleep on the floor, 1 step | 1,689.4 | 1,241.2 | 1.36 |
| step: 60 chains of 8 links on distance, ball and hinge joints, 30 steps | 208,366.5 | 175,390.4 | 1.19 |
| step: 200 swept spheres fired at walls, spheres, compounds and terrain, swept and pruned, 6 steps | 131,362.8 | 124,125.4 | 1.06 |
| step: 200 swept spheres fired at walls, spheres and compounds in the grid, 6 steps | 71,724.5 | 61,023.7 | 1.18 |
| query: 2,048 rays in a batch, as columns | 58,853.1 | 39,985.3 | 1.47 |
| query: 2,048 rays in a batch, as columns, unbounded and downward | 45,144.6 | 32,774.6 | 1.38 |
| query: 2,048 rays in a batch, as records | 56,696.8 | 37,858.8 | 1.50 |
| query: 512 rays one by one | 14,572.8 | 9,616.3 | 1.52 |
| query: 1,024 sphere casts in a batch, as columns | 19,071.9 | 13,541.2 | 1.41 |
| query: 2,048 boxes of a few cells | 28,116.1 | 17,824.2 | 1.58 |
| query: 64 boxes of more cells than the grid is walked for | 39,415.8 | 22,473.7 | 1.75 |
| drain: transforms of 1,001 bodies | 60.8 | 40.9 | 1.48 |
| drain: linear state of 1,000 dynamic bodies, with user data and without handles | 74.2 | 61.6 | 1.20 |
| drain: contacts of the pile, as columns | 94.8 | 90.1 | 1.05 |
| drain: contacts of the pile, as records | 193.9 | 191.9 | 1.01 |
| drain: contact events of the pile, without the ones that stay | 66.0 | 63.5 | 1.04 |
| drain: stats of the pile | 33.6 | 19.0 | 1.77 |
| build: a mesh of 8,192 triangles from tables | 2,835.6 | 2,702.9 | 1.05 |
| build: a mesh of 8,192 triangles from buffers | 2,346.8 | 2,231.3 | 1.05 |
| build: a mesh of 1,408 triangles without indices | 612.7 | 548.7 | 1.12 |
| build: a heightfield 65 by 65 from a table | 55.8 | 39.2 | 1.42 |
| build: a heightfield 65 by 65 from a buffer | 18.1 | 10.0 | 1.80 |
| build: add 5,000 bodies of every primitive, then remove every other one | 9,485.2 | 8,946.3 | 1.06 |
| sync: move the origin under 5,000 bodies, and back | 832.6 | 617.1 | 1.35 |
| sync: apply the transforms of 5,000 bodies | 2,488.0 | 2,296.4 | 1.08 |

**Read the floor first.** The rows nothing touched (contacts drained as columns and as records, contact events, adding and removing bodies, applying transforms) read 1.01 to 1.08 in this run, and 0.98 to 1.03 in the same A/B run half an hour earlier on the commit before the last. That is the noise: four other agents' tests were running on the machine. A row inside it is "no change": those, the swept spheres against terrain (1.06 here, 1.01 before), and the three mesh builds (1.05, 1.05 and 1.12 here; 1.00, 0.98 and 1.02 before).

**What moved, and why.** In the order the commits were made; the ratios in passing are from three-round A/Bs taken as each landed.

- *Piles and sweeps, 1.09x to 1.38x at once: the pair filter (`b1e983d`).* A broadphase offers a pair, and most pairs it offers are apart: 21,000 offers a step for the 1,000 spheres, of which 3,500 go on to the narrow phase. 0.3.2 checked liveness, a joint table, body types and masks, and read both shapes, before it compared any bounds. It now compares the bounds first, from a column of sphere radii the world keeps beside its shapes, and reads nothing else for a pair that is apart. The pairs that survive are the same ones in the same order.
- *Piles, another tenth or so (`071ce44`).* Contact preparation and the solver read forty-odd columns from the world for every contact in every iteration. They now bind each column once, with record patterns, and the cache match takes the anchors preparation has just made in place of making them again.
- *The 5,000 drifting spheres, 1.03x to 1.18x (`34c32e9`), then to 1.3x (`ba063aa`).* A body crossing a cell rebuilds the grid's occupied-neighbour pairs: thirteen lookups a cell, nearly all of them misses. A byte for each residue of a packed cell now says a cell is empty without the lookup. The rest is the same binding of columns in the integrator, the island builder, the sleep passes, the contact cache and the drains, which is also where the transform and linear-state drains (1.48x, 1.2x), the origin shift (1.35x) and the sleeping floor (1.36x) come from.
- *Rays, sphere casts and box queries, 1.38x to 1.75x (`52e2bb1`).* A ray walked the grid by probing all twenty-seven cells around each cell it entered, eighteen of which it had probed from the cell before. It now probes the whole neighbourhood once and after that only the cells a step brings in; a body is still tested once, at the same point in the walk, against the same best distance. Box queries pack a row's cell keys by adding one, ask the filter before a lookup, and drop a set of marks that could never fire, a body being in one cell.
- *A sweep that starts disordered, 6.4x (`640ec2e`).* The sweep-and-prune order is insertion-sorted, which is right when motion is coherent and quadratic when 2,000 bodies have just been added. Once the sort has moved more than eight times the bodies there are, what is left is handed to `data.argsort` over a buffer of keys, filled by a buffer sink and applied by a comprehension over a typed slice of the permutation. Both sorts are stable, so they end in the same order; keys that are not numbers, which the two would place differently, keep the insertion sort.
- *Compounds, 1.17x to 1.26x (`ba0c230`).* Two three-piece compounds are nine piece pairs, and each pair placed both its pieces in the world again. A side's pieces are now placed once.
- *Heightfields, 1.4x from a table and 1.8x from a buffer (`c32b927`, `819f278`).* The least and greatest sample are one `min` and one `max` on a typed receiver. That took two tries: see below.
- *`stats`, 1.77x (`99d9002`).* Two `#[...]` counts over columns bound once, where two loops read the world for every body and contact.

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
| step: 600 spheres, capsules and boxes on a heightfield 65 by 65, 30 steps, which 0.3.2 gets wrong | 177,894.5 |
| step: 600 spheres, capsules and boxes on a mesh of 8,192 triangles, 30 steps, which 0.3.2 gets wrong | 296,480.4 |
| step: 200 swept spheres fired at walls, spheres, compounds and terrain, 6 steps, which 0.3.2 gets wrong | 124,222.1 |

The swept spheres compute, under the grid, exactly what they compute swept and pruned: the same contacts and the same checksum.

A fourth defect was fixed on the way (`c62e8c7`): a body refused for its rotation had already taken a slot, and left a live body with no rotation behind for the next step to trip over.
