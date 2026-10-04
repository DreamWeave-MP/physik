# physik.jpk

`physik` is a pure-Luau 3D rigid-body kernel for Jess. It uses Rapier and Jolt as architectural references, but the implementation is shaped around Jess/Luau's actual cost model instead of transliterating Rust or C++.

The package exports:

```luau
local physik = require("@dream/physik")
```

## Design shape

The hot path is deliberately data-oriented:

- bodies are generation-checked numeric handles, not object tables;
- hot translational/rotational state lives in parallel arrays of Luau native `vector` values;
- scalar metadata lives in parallel arrays and is never materialized as a body object during stepping;
- contacts live in reusable parallel arrays;
- body rotation axes are cached once per collision phase and reused by OBB SAT and inverse inertia;
- sphere inertia has an isotropic fast path;
- broad phase defaults to a numeric 3D hash grid for cell-sized bodies plus a coherent X sweep for bodies too large for one cell;
- a pure sweep-and-prune broad phase remains selectable for comparison/special workloads;
- contact normal and two fixed friction tangent effective masses are prepared once before velocity iterations.

No per-contact table is allocated while stepping. `world:getBody`, `world:contacts`, and ray-query result tables are intentionally boundary APIs, not the simulation representation.

## Shapes

```luau
physik.sphere(radius)
physik.box(vector.create(halfX, halfY, halfZ))
physik.capsule(radius, halfHeight) -- local Y axis
```

Rigid-body collision coverage in v0.1:

- sphere / sphere
- sphere / oriented box
- oriented box / oriented box using 15-axis SAT
- capsule / sphere
- capsule / capsule
- capsule / box

Capsule/box currently probes the two capsule endpoints plus midpoint. This is deliberately compact rather than a complete segment-vs-OBB manifold implementation.

## Minimal example

```luau
local physik = require("@dream/physik")

local world = physik.newWorld({
    gravity = vector.create(0, -9.81, 0),
    velocityIterations = 8,
})

world:addBody({
    type = "static",
    shape = physik.box(vector.create(50, 0.5, 50)),
    position = vector.create(0, -0.5, 0),
})

local body = world:addBody({
    shape = physik.sphere(0.5),
    position = vector.create(0, 10, 0),
    mass = 1,
    friction = 0.6,
    restitution = 0.05,
})

for _ = 1, 600 do
    world:step(1 / 60)
end

print(world:getBody(body).position)
```

## Body API

`world:addBody(desc)` accepts:

```text
type             dynamic (default), static, or kinematic
shape            sphere(0.5) by default
position         vector
rotation         physik.quat(...) or physik.axisAngle(...)
velocity         vector
angularVelocity  vector
mass             dynamic bodies; default 1
friction         default 0.6
restitution      default 0.05
linearDamping    default 0.02
angularDamping   default 0.05
group            32-bit collision group bit(s), default 1
mask             32-bit collision mask, default all bits
userData         arbitrary caller value
```

Mutation/query surface:

```text
world:removeBody(handle)
world:getBody(handle)
world:setPosition(handle, position [, wake])
world:setRotation(handle, quat [, wake])
world:setVelocity(handle, velocity)
world:setAngularVelocity(handle, velocity)
world:wake(handle)
world:applyForce(handle, force)
world:applyForceAtPoint(handle, force, worldPoint)
world:applyImpulse(handle, impulse [, worldPoint])
```

Handles encode slot + generation. Removing and reusing a slot does not make an old handle valid again.

## Constraints

v0.1 includes a compact distance constraint:

```luau
local joint = world:addDistanceJoint(a, b, {
    localAnchorA = vector.create(0, 0, 0),
    localAnchorB = vector.create(0, 0, 0),
    length = 2,
    stiffness = 0.25,
    damping = 0.05,
})

world:removeJoint(joint)
```

## Queries

```luau
local hit = world:raycast(origin, direction, maxDistance, mask)
```

The result contains `handle`, `distance`, `point`, `normal`, and `userData`.

`world:contacts()` materializes the current contacts only when requested. `world:stats()` reports frame, body, sleeping, candidate, contact, and joint counts.

## World tuning

```text
gravity                 default (0, -9.81, 0)
velocityIterations      default 8
positionCorrection      default 0.65
penetrationSlop         default 0.003
restitutionThreshold    default 1.0
allowSleep              default true
sleepLinear             default 0.04
sleepAngular            default 0.05
sleepDelay              default 0.6
wakeSpeed               default 0.5
wakePenetration         default 0.05
broadphase              "grid" (default) or "sap"
broadphaseCellSize      default 1.0
```

Sleeping bodies are solver-static until meaningful motion or penetration wakes them. This avoids routine support impulses constantly waking settled bodies. v0.1 does not yet build explicit simulation islands, so very large mutually-supported piles can take longer to become fully quiescent than a production island sleeper.

## Measured sandbox performance

These are implementation-development measurements from the supplied Jess binary in the ChatGPT Linux sandbox, not universal promises. They are useful for regression direction only.

Unprofiled, 60 Hz, one substep:

```text
1000 dynamic spheres + static floor, 120 steps:
~0.660 s total
~5.50 ms / step
~182k body-steps / second

250 dynamic oriented boxes + static floor, 120 steps:
~0.538 s total
~4.48 ms / step
~55.8k body-steps / second
```

The box scene becomes substantially more expensive once dense OBB contacts accumulate. Jess's sampling/native-block profiler adds material overhead to this workload, so raw and profiled timings should not be compared directly.

Important optimizations discovered by profiling:

- numeric 3D grid beats single-axis SAP for compact piles;
- a one-cell default near the common collider diameter beat the original 1.25-cell default in the sphere pile;
- oversized/rotated bodies must not fall into an O(L²) fallback; they now use the coherent sweep path;
- allocating a nested SAT axis-test closure per box/box candidate caused tens of MB of transient growth; flattening it removed almost all of that growth;
- caching box axes avoids recomputing quaternion bases per candidate;
- isotropic sphere inverse inertia bypasses all orientation transforms;
- fixed contact tangent bases let both friction effective masses be hoisted out of solver iterations.

## Deliberate v0.1 limits

This is already a real rigid-body kernel, but it is not pretending to be Jolt or Rapier feature-complete.

Not implemented yet:

- persistent 2-4 point contact manifolds;
- cross-frame contact impulse warm starting;
- explicit simulation-island construction/sleeping;
- CCD / TOI solving;
- convex hulls, triangle meshes, heightfields, compound shapes;
- hinge/slider/cone/motor constraints;
- exact capsule-vs-box manifold generation;
- broad-phase acceleration for raycasts;
- deterministic fixed-point simulation.

The next quality jump should be **persistent manifolds + warm starting**, followed by **islands**, before adding a pile of exotic shape types.
