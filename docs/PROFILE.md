# physik profiling / L3i pressure notes

This file records where the pure-Luau implementation actually spent time while being developed against the supplied Jess binary.

## What Luau is already surprisingly good at

Native `vector` values are excellent hot-state storage. A direct array of vectors beat reconstructing vectors from packed `buffer` f32 triples in the initial microbenchmark. That is why positions, velocities, forces, angular velocities, quaternion vector parts, AABBs, and cached axes remain vector arrays instead of being packed blindly.

The useful split is not "pack everything". It is:

```text
hot math state      native vector / number arrays
cold identity/meta  compact scalar state where profitable
boundary objects    allocate only when the caller asks
```

## Hot spots after the current optimization pass

For dense oriented-box scenes the remaining time is primarily:

1. 15-axis OBB SAT;
2. broad-phase pair enumeration / AABB rejection;
3. contact impulse application and inverse inertia;
4. preparation of three effective masses per contact.

For sphere piles, after the isotropic-inertia fast path, broad-phase enumeration is the dominant cost.

## Native primitives that would move the needle

These are candidates for L3i only if they can stay generic rather than becoming `physik`-specific host policy.

### 1. Typed native arenas / structure-of-arrays

A generic growable typed arena for numbers, vectors, u32/i32, and maybe packed records would let JPKs keep stable integer handles while avoiding Lua table slot overhead.

Desirable shape:

```text
arena.vector3(capacity)
arena.f32(capacity)
arena.u32(capacity)
get/set by integer index
cheap grow/reserve
no userdata allocation per element
```

This would help physics, ECS, scenegraphs, mesh processing, navigation, and SQL-like column stores at once.

### 2. Batched spatial hash primitive

A generic numeric-key bucket build from `(key, value)` integer arrays could erase much of the per-frame table work in the broad phase.

Useful operations:

```text
clear/rebuild bucket heads
append integer value to numeric key
iterate bucket
```

No physics semantics required.

### 3. Quaternion / basis transforms

Native helpers for:

```text
quat_rotate(q, v)
quat_inverse_rotate(q, v)
quat_to_basis(q)
diag_basis_transform(basis, diagonal, v)
```

would hit physics and scenegraph/animation workloads. The current implementation already caches bases aggressively because recomputing them is visible in profiles.

### 4. Small fixed vector batches

A primitive that evaluates several dot/cross projections without table creation could help SAT and geometry kernels, but only add this if other JPKs show the same pressure. A physics-only "obb_sat" intrinsic would be the wrong layer.

### 5. Sort / radix-sort over numeric side arrays

A stable native sort for numeric keys with integer payloads would improve sweep broad phases, render sorting, ECS batching, and database/order-by workloads.

Candidate:

```text
sort_pairs(keys, values, count)
```

or a reusable permutation output.

## What *not* to native-ize yet

Do not add a native rigid-body engine API merely because this JPK exists. The point of this package is to establish how far typed Luau can go and to expose generic substrate weaknesses. Persistent manifolds, islands, CCD, and more constraints should be implemented/prototyped in Luau first unless profiling demonstrates a reusable missing primitive.
