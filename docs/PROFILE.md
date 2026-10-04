# physik profiling / L3i pressure notes

This records what the pure-Luau 3D physics implementation actually revealed while being developed against the supplied Jess binary.

## What Luau already does well

Native `vector` values are excellent hot-state storage. A direct vector array beat reconstructing vectors from packed `buffer` f32 triples in the initial layout microbenchmark.

The useful split remains:

```text
hot math state      native vector / number arrays
cold identity/meta  compact scalar state where profitable
boundary objects    allocate only when the caller asks
```

Blindly packing hot vectors is slower here because decode/reconstruction cost dominates the memory win.

## Important profiling discoveries

### Hybrid broad phase needs an adaptive large-body path

The numeric 3D grid is very good for compact bodies, but rotated boxes frequently have AABBs larger than one cell. An early fallback performed O(L^2) work over oversized bodies; replacing that with a coherent X sweep removed that cliff.

A second pathology appeared later: **one enormous floor collider** enabled the all-body sweep, causing it to scan many small-small pairs that the grid had already handled. v0.2 therefore chooses:

```text
few oversized bodies   direct large-vs-world checks
many oversized bodies  coherent X sweep
```

That change moved the 1000-sphere + large-floor benchmark back below the original v0.1 cost even after adding manifolds, cache matching, and islands.

### Inner-loop closures are poison

The first OBB SAT implementation allocated a nested axis-test closure per box/box candidate. In a 250-box tumble, profiler-observed memory growth was roughly 31 MB. Flattening that helper removed nearly all of that transient growth (about 1.7 MB in the comparable sampled run) and materially improved runtime.

### Precompute angular response, not only scalar effective mass

For a contact axis `n`, the solver repeatedly needs:

```text
I^-1 (r x n)
```

v0.2 prepares this angular response for the normal and both friction tangents once per active contact, then reuses the vectors directly inside every solver iteration. Recomputing quaternion/basis inertia transforms inside each impulse application is gratuitous in Luau.

### Warm starting changes the stabilization regime

Once 2-4 point manifolds and cross-frame impulse matching exist, the old v0.1 defaults were too aggressive. Large Baumgarte velocity bias plus explicit position correction injects visible energy into dense multi-point support contacts.

The current default is therefore:

```text
velocityIterations = 4
baumgarte = 0.03
positionCorrection = 0.12
warmStarting = true
warmStartFactor = 0.85
```

In the 250-box / 300-step regression scene, four warm-started iterations showed lower residual motion than eight cold iterations while also running faster.

### Sleeping must remove work, not just zero velocity

The first island pass correctly slept connected bodies but still prepared every sleeping contact and walked it in every solver iteration. v0.2 retains those contacts for wake detection while marking them inactive for preparation/solve.

1000 independent resting spheres improved from about 3.62 ms/step to roughly 1.84-2.15 ms/step in repeated sandbox runs after this cut.

### Exact capsule/box was worth doing in Luau

The original capsule/box narrow phase sampled only the two capsule endpoints plus its midpoint. That is cheap but misses a small box touching a long capsule between those samples. v0.2 now transforms the capsule segment into box-local space, tests segment/AABB intersection, and otherwise finds the closest feature against all twelve box edges without allocation. A regression test places a box one unit along a long capsule spine specifically where the old three-sample path would miss it.

### Sphere-first CCD is a useful correctness slice

Opt-in dynamic spheres now retain their previous position and sweep against static/kinematic spheres and boxes before the discrete pass. This is deliberately not a general TOI engine, but it closes the most embarrassing tunneling case: a fast projectile crossing a thin wall in one 60 Hz tick. Because the path is opt-in and scans targets only for marked spheres, it does not tax the ordinary pile benchmarks.

### Connected constraints need collision semantics

Ball/socket and hinge constraints initially exposed a misleading failure mode when the two connected shapes overlapped: the constraint solver and contact solver fought each other. v0.2 now maintains a numeric body-pair suppression count for joints and defaults `collideConnected` to false. Multiple joints between one pair are reference-counted, so removing one constraint does not accidentally re-enable contact while another still exists.

This is small policy, but it matters for correctness far more than adding another exotic joint type.

### Tiny-table shape matters more than slogans

The persistent manifold cache went through three representations in the same 250-box profiled workload:

```text
nested pair entry + five child arrays     ~3.38 MB bytes grown
global SoA columns keyed by cache slot    ~4.09 MB bytes grown
one pre-sized numeric array per pair      ~2.96 MB bytes grown
```

The winning entry is a single `table.create(26)` addressed entirely by integer offsets: six metadata slots plus four fixed-width contact records. It removes five child-table allocations per persistent pair without forcing many independently-growing global arrays.

This is a useful warning for other JPKs: **measure the VM's actual table-growth/allocation behavior before treating SoA as automatically more compact.** Typed native arenas would likely change this tradeoff again.

## Current hot spots

For dense oriented-box scenes:

1. 15-axis OBB SAT and face clipping;
2. broad-phase pair enumeration / AABB rejection;
3. active-contact preparation;
4. sequential normal/friction iterations.

For compact sphere piles, broad-phase enumeration remains dominant once inertia transforms disappear.

## Generic L3i primitives that would move the needle

These are useful only if they remain generic substrate rather than becoming physics host policy.

### 1. Typed native arenas / SoA columns

A growable typed arena for numbers, vectors, u32/i32, and perhaps fixed-width records would let JPKs keep stable integer handles without Lua table-slot overhead.

Desirable shape:

```text
arena.vector3(capacity)
arena.f32(capacity)
arena.u32(capacity)
indexed get/set
reserve/grow
no userdata allocation per element
```

This would help physics, ECS, scenegraphs, mesh processing, navigation, and SQL-like column stores simultaneously.

### 2. Numeric spatial-hash buckets

A reusable native numeric-key bucket builder would remove much of the per-frame table work from the broad phase:

```text
clear/rebuild bucket heads
append integer payload to numeric key
iterate bucket contents
```

No physics semantics required.

### 3. Quaternion / basis transforms

Generic helpers for:

```text
quat_rotate(q, v)
quat_inverse_rotate(q, v)
quat_to_basis(q)
diag_basis_transform(basis, diagonal, v)
```

would help physics and scenegraph/animation workloads. Physik already caches bases aggressively because recomputation is visible.

### 4. Stable numeric sort / radix permutation

A native stable sort over numeric keys with integer payloads would improve sweep broad phases, render sorting, ECS batching, and database/order-by workloads.

Potential shape:

```text
sort_pairs(keys, values, count)
```

or a reusable permutation output.

### 5. Compact union-find arena

Island construction is already cheap enough in Luau for current loads, but a generic integer union-find over typed arrays could be useful to ECS grouping, graph connectivity, and physics alike if those workloads converge on the same primitive.

## What not to native-ize yet

Do **not** add `physics_step`, `obb_sat`, or a native rigid-body engine API merely because Physik exists. The point of this package is to identify generic substrate weaknesses while keeping simulation policy in Luau.

CCD, meshes, additional constraints, and island scheduling should continue in Luau first. Only promote the repeated, reusable pressure points into L3i.
