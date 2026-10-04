# physik.jpk

`physik` is a pure-Luau 3D rigid-body kernel for Jess. Rapier and Jolt are architectural references; the implementation is intentionally shaped around Jess/Luau's actual cost model instead of transliterating Rust or C++.

The package exports:

```luau
local physik = require("@dream/physik")
```

Current package version: **0.2.0**.

## v0.2 architecture

The hot path is data-oriented and allocation-averse:

- generation-checked numeric body handles rather than per-body simulation objects;
- native Luau `vector` arrays for hot translational and rotational state;
- parallel scalar arrays for metadata;
- reusable SoA contact storage;
- cached per-frame body bases reused by OBB SAT and inverse inertia;
- isotropic sphere inertia fast path;
- numeric 3D hash-grid broad phase for cell-sized bodies;
- adaptive oversized-body handling: a few large bodies use direct large-vs-world checks, while dense oversized sets use a coherent X sweep;
- 2-4 point clipped box-face manifolds;
- cross-frame contact matching in local body space;
- persistent-pair cache entries packed into one pre-sized integer-indexed table (no child contact tables);
- cached accumulated normal/friction impulses and warm starting;
- fixed contact tangent bases and precomputed normal/tangent effective masses;
- precomputed angular response vectors reused directly by the iterative solver;
- union-find simulation islands spanning dynamic contact/joint graphs;
- island sleeping/waking;
- sleeping contacts are retained for wake detection but removed from contact preparation and solving.

No per-contact table is allocated while stepping. `world:getBody`, `world:contacts`, ray-query results, and diagnostics intentionally materialize tables only at API boundaries.

## Shapes

```luau
physik.sphere(radius)
physik.box(vector.create(halfX, halfY, halfZ))
physik.capsule(radius, halfHeight) -- local Y axis
```

Collision coverage:

- sphere / sphere
- sphere / oriented box
- oriented box / oriented box using 15-axis SAT
- capsule / sphere
- capsule / capsule
- capsule / box

Box/box face contacts are clipped to a manifold of at most four spatially spread points. Edge/edge SAT contacts remain single-point. Capsule/box uses a full segment-vs-OBB closest-feature query, including box-edge tests; it is no longer limited to endpoint/midpoint probes.

## Minimal example

```luau
local physik = require("@dream/physik")

local world = physik.newWorld({
    gravity = vector.create(0, -9.81, 0),
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
ccd              opt-in swept CCD for dynamic spheres; default false
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

Handles encode slot + generation. Removing and reusing a slot never makes an old handle valid again.

## Constraints

v0.2 includes distance, ball/socket, and hinge constraints. All participate in dynamic-island construction. Connected-body collision is disabled by default; set `collideConnected = true` on the joint descriptor to keep it enabled.

Distance:

```luau
local joint = world:addDistanceJoint(a, b, {
    localAnchorA = vector.create(0, 0, 0),
    localAnchorB = vector.create(0, 0, 0),
    length = 2,
    stiffness = 0.25,
    damping = 0.05,
})
```

Ball/socket:

```luau
local joint = world:addBallJoint(a, b, {
    localAnchorA = vector.create(0, 0, 0),
    localAnchorB = vector.create(0, 1, 0),
    stiffness = 0.35,
    damping = 0.1,
})
```

Hinge with an optional torque-limited motor:

```luau
local joint = world:addHingeJoint(a, b, {
    localAnchorA = vector.create(0, 0, 0),
    localAnchorB = vector.create(0, 1, 0),
    localAxisA = vector.create(0, 1, 0),
    localAxisB = vector.create(0, 1, 0),
    angularStiffness = 0.25,
    motorSpeed = 3,
    maxMotorTorque = 20,
})
```

Remove any joint with:

```luau
world:removeJoint(joint)
```

`collideConnected` defaults to false for all three joint types. The world maintains an integer pair-suppression count so multiple joints between the same body pair compose correctly. Joint removal, stale-joint cleanup, and body removal all release their suppression entries.

## Queries and diagnostics

```luau
local hit = world:raycast(origin, direction, maxDistance, mask)
```

The result contains `handle`, `distance`, `point`, `normal`, and `userData`.

`world:contacts()` materializes current contacts on demand. `world:stats()` currently reports:

```text
frame
bodies
sleeping
islands
candidates
contacts
activeContacts
manifolds
warmStarted
ccdHits
joints
```

`contacts` includes sleeping support contacts retained for wake detection. `activeContacts` is the subset actually prepared/solved this step.

## World tuning

```text
gravity                    default (0, -9.81, 0)
velocityIterations         default 4
positionIterations         default 1
positionCorrection         default 0.12
baumgarte                  default 0.03
baumgarteSlop              default 0.05
maxDepenetrationVelocity   default 1.0
penetrationSlop            default 0.003
restitutionThreshold       default 1.0
allowSleep                 default true
sleepLinear                default 0.04
sleepAngular               default 0.05
sleepDelay                 default 0.6
wakeSpeed                  default 0.5
wakePenetration            default 0.05
warmStarting               default true
warmStartFactor            default 0.85
manifoldMatchDistance      default 0.08
contactMargin              default 0.015
broadphase                 "grid" (default) or "sap"
broadphaseCellSize         default 1.0
directLargeThreshold       default 16
ccdSlop                    default 0.001
```

The v0.2 defaults changed deliberately. Multi-point manifolds plus cross-frame warm starting made the v0.1 combination of eight iterations and aggressive Baumgarte velocity bias unnecessarily energetic. v0.2 lowers the default to four iterations and a 0.03 Baumgarte factor while retaining explicit positional correction.

## Continuous collision detection

v0.2 includes an intentionally narrow first CCD path. Set `ccd = true` on a **dynamic sphere** body to sweep its center from the previous position to the integrated position before the discrete collision pass.

Current sweep targets:

- static or kinematic spheres;
- static or kinematic oriented boxes.

The sphere/box sweep uses the box expanded by the moving sphere radius, making it conservative near box corners but effective at preventing tunneling through thin world geometry. On a hit Physik rewinds to just before TOI using `ccdSlop`, removes the inward normal velocity with restitution, and then lets the ordinary discrete solver continue.

Current CCD does **not** handle dynamic-vs-dynamic sweeps, capsules, boxes as moving sweep shapes, or rotational TOI. It is opt-in so ordinary bodies do not pay the O(N) target scan. `world:stats().ccdHits` reports the number of CCD clamps in the most recent substep.

## Measured sandbox performance

These are development measurements from the supplied Jess binary in the ChatGPT Linux sandbox. They are regression-direction measurements, not hardware-independent promises.

Unprofiled, 60 Hz, one substep, v0.2 defaults / four velocity iterations:

```text
1000 dynamic spheres + large static floor, 120 steps:
4.62 / 4.67 / 4.65 ms per step across three runs
~215k body-steps / second typical

250 dynamic oriented boxes + static floor, 120 steps:
5.00 / 4.98 / 4.79 ms per step across three runs
~50-52k body-steps / second

1000 independent spheres fully asleep on a floor, 300 measured steps:
2.15 / 1.95 / 1.84 ms per step across three runs
```

The sleeping-world optimization is significant because the engine still generates the 1000 support contacts each frame for wake semantics, but none are prepared or iterated once their islands are asleep. Before that cut, the same settled scene measured about **3.62 ms/step** in this sandbox.

A 250-box, 300-step convergence comparison with the v0.2 stabilization settings produced:

```text
4 iterations + warm starting:
~6.55 ms/step, residual-motion proxy 17.9

8 iterations + warm starting disabled:
~6.95 ms/step, residual-motion proxy 70.9
```

That benchmark is intentionally chaotic and the residual-motion value is only a regression metric, but it demonstrates the intended direction: persistent impulses can buy convergence rather than merely adding bookkeeping.

## Deliberate v0.2 limits

This is now a real compact rigid-body kernel, not a claim of Rapier/Jolt feature parity.

Not implemented yet:

- general CCD / TOI (dynamic-vs-dynamic, moving boxes/capsules, rotational sweeps);
- convex hulls, triangle meshes, heightfields, compound shapes;
- slider / cone / six-DOF constraints;
- edge/edge multi-point manifolds;
- broad-phase acceleration for raycasts;
- deterministic fixed-point simulation;
- solver parallelism across independent islands.

The next major correctness step should extend CCD beyond sphere-vs-static/kinematic targets. The next major shape feature should probably be **triangle meshes / heightfields**, because those matter far more to game-world usefulness than exotic convex primitives.
