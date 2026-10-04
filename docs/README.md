# physik.jpk

`physik` is a pure-Luau 3D rigid-body kernel for Jess. Rapier and Jolt are architectural references; the implementation is deliberately shaped around Jess/Luau's actual cost model rather than transliterating Rust or C++.

```luau
local physik = require("@dream/physik")
```

Current package version: **0.3.0**.

## v0.3 architecture

The hot path is data-oriented and allocation-averse:

- generation-checked numeric handles instead of per-body simulation objects;
- native Luau `vector` arrays for hot position/velocity/angular state;
- parallel scalar arrays for metadata;
- reusable SoA contact storage;
- cached body bases reused by SAT and inverse-inertia work;
- isotropic sphere inertia fast path;
- numeric 3D hash-grid broadphase for compact bodies;
- adaptive oversized-body path: direct large-vs-world when only a few large bodies exist, coherent X sweep when many do;
- 2-4 point clipped box-face manifolds;
- local-anchor contact matching across frames;
- fixed-width integer-indexed persistent contact-cache records;
- accumulated normal + 2D friction warm starting;
- precomputed normal/tangent effective masses and angular response vectors;
- union-find simulation islands spanning contacts and joints;
- island-wide sleep/wake;
- sleeping contacts retained for wake detection but removed from active solve work;
- buffer-backed triangle meshes and heightfields;
- dynamic compound colliders assembled from spheres, boxes, and capsules;
- lazy spatial-query indexing shared with the broadphase grid;
- reusable SoA drains for transforms, contacts, contact-pair events, ray batches, and sphere-cast batches.

Ordinary stepping does not allocate per-contact tables. Materialized records are limited to convenience APIs such as `getBody`, `contacts`, single `raycast`, and single `sphereCast`.

## Shapes

```luau
physik.sphere(radius)
physik.box(vector.create(halfX, halfY, halfZ))
physik.capsule(radius, halfHeight) -- local Y axis
physik.compound(children)
physik.triangleMesh(vertices, indices [, opts])
physik.heightfield(width, depth, heights [, opts])
```

`triangleMesh` accepts either vector/u32 tables or packed buffers. Vertex buffers are packed f32 xyz triples. Index buffers are packed u32 triangle lists and are zero-based when supplied as a buffer. Table indices are one-based.

Triangle meshes and heightfields are static/kinematic world geometry; Physik intentionally rejects dynamic triangle soup. Use compound colliders for dynamic ships/actors made of several primitive pieces.

Primitive collision coverage includes sphere/sphere, sphere/OBB, OBB/OBB 15-axis SAT, capsule/sphere, capsule/capsule, exact capsule/OBB closest-feature tests, primitive/triangle tests, and compound traversal. World-geometry contacts are capped to four deepest distinct points per body pair.

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
})

for _ = 1, 600 do
    world:step(1 / 60)
end

print(world:getBody(body).position)
```

## Bodies and ECS integration

`world:addBody(desc)` supports dynamic (default), static, and kinematic bodies. Important descriptor fields:

```text
type
shape
position
rotation
velocity
angularVelocity
mass
friction
restitution
linearDamping
angularDamping
group
mask
userData
sensor
ccd
```

Core mutation surface:

```text
world:removeBody(handle)
world:setPosition(handle, position [, wake])
world:setRotation(handle, quat [, wake])
world:setVelocity(handle, velocity)
world:setAngularVelocity(handle, velocity)
world:driveKinematic(handle, targetPosition [, targetRotation], dt)
world:wake(handle)
world:applyForce(handle, force)
world:applyForceAtPoint(handle, force, worldPoint)
world:applyTorque(handle, torque)
world:applyImpulse(handle, impulse [, worldPoint])
world:applyAngularImpulse(handle, impulse)
```

Bulk ECS synchronization avoids one record per body:

```luau
local state, count = world:writeTransforms(reusable, {
    awakeOnly = true,
    dynamicOnly = true,
})

-- state.handles
-- state.positions
-- state.rotations       -- quaternion xyz as vector
-- state.rotationW
-- state.velocities
-- state.angularVelocities

world:applyTransforms(state) -- authoritative/replicated write-back
```

`applyTransforms` invalidates persistent contact anchors and the lazy spatial-query index, so an authoritative teleport is visible to collision queries immediately.

For large-coordinate worlds, `world:shiftOrigin(delta)` rebases every live body and current contact point. The numeric hash grid is intentionally optimized for a local simulation bubble rather than pretending arbitrary astronomical coordinates are free.

## Fixed stepping

`world:step(dt [, substeps])` performs one explicit simulation step.

`world:advance(realDt [, fixedDt [, maxCatchUpSteps]])` maintains an accumulator and returns:

```text
stepsExecuted, interpolationAlpha, droppedCatchUpTime
```

This is the intended host-facing path when render/network time and simulation tick time differ.

## Constraints

v0.3 provides:

- distance joints;
- ball/socket joints;
- hinge joints;
- optional torque-limited hinge motors.

```luau
world:addDistanceJoint(a, b, opts)
world:addBallJoint(a, b, opts)
world:addHingeJoint(a, b, opts)
world:removeJoint(id)
```

Connected-body collision is disabled by default through reference-counted numeric pair suppression. Set `collideConnected = true` on a joint when contact between the connected bodies is intentional.

## Queries

Single ray:

```luau
local hit = world:raycast(origin, direction, maxDistance, mask)
```

Reusable batches:

```luau
world:raycastBatch(origins, directions, maxDistances, masks, reusable)
world:raycastBatchSoA(origins, directions, maxDistances, masks, reusable)
```

Spatial overlap:

```luau
local handles, count = world:queryAabb(minimum, maximum, mask, reusable)
```

Sphere sweep / shape cast:

```luau
local hit = world:sphereCast(startPos, endPos, radius, mask, ignoreHandle)
world:sphereCastBatchSoA(starts, ends, radii, masks, ignoreHandles, reusable)
```

Raycasts and AABB queries lazily reuse Physik's body grid. Large bodies are checked directly, compact bodies are traversed spatially, and pathological huge query ranges fall back to a linear scan rather than walking unbounded cell counts.

Sphere casts reuse the same TOI kernels as CCD against spheres, boxes, capsules, compounds, triangle meshes, and heightfields.

## Contacts, sensors, and gameplay events

Current contact state can be drained without per-contact objects:

```luau
local contacts, count = world:writeContacts(reusable)
```

For triggers/gameplay pair transitions, enable event tracking explicitly:

```luau
local world = physik.newWorld({ contactEvents = true })
...
local events, count = world:writeContactEvents(reusable)
```

`events.kind[n]` is one of:

```luau
physik.contactEvent.enter
physik.contactEvent.stay
physik.contactEvent.exit
```

`writeContactEvents(out, false)` omits `stay` records. Event tracking is **off by default** so worlds that do not need gameplay pair transitions do not pay that bookkeeping cost.

Sensor bodies participate in contact/event detection but never enter islands or impulse solving.

## Continuous collision detection

Set `ccd = true` on a **dynamic sphere** to sweep it from its previous position to its integrated position before the discrete contact pass.

Current CCD targets:

- static/kinematic spheres;
- oriented boxes;
- capsules;
- compound primitive children;
- triangle meshes, including face/edge/vertex TOI;
- heightfields.

CCD remains intentionally sphere-first. There is no dynamic-vs-dynamic TOI island solver, moving-box/capsule cast, or rotational TOI yet. Unsupported fast movers should use substeps or explicit sphere casts where appropriate.

## World tuning

Important options include:

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
ccdSlop                    default 0.001
broadphase                 "grid" (default) or "sap"
broadphaseCellSize         default 1.0
directLargeThreshold       default 16
meshQueryCellLimit         default 4096
queryRayCellLimit          default 2048
queryAabbCellLimit         default 8192
fixedDt                    default 1/60
maxCatchUpSteps            default 8
contactEvents              default false
```

## Measured sandbox performance

These are regression-direction samples from the supplied Jess binary in the ChatGPT Linux sandbox, not hardware-independent promises. Recent unprofiled v0.3 samples:

```text
1000 dynamic spheres + static floor, 120 steps:
~5.37 ms/step

250 tumbling OBBs + floor, 120 steps:
~5.53 ms/step

256 three-piece dynamic compound "ships" + floor, 120 steps:
~4.42 ms/step

65x65 heightfield + 400 active spheres:
~3.12 ms/step

same terrain as an 8192-triangle generic mesh + 400 active spheres:
~4.60 ms/step
```

The body-query grid matters more as scene population grows. In a 4,096-small-body scene with one large floor, a 20,480-ray batch sample measured approximately:

```text
forced linear body scan:  ~372 us/ray
spatial grid traversal:   ~12.5 us/ray
```

That is roughly a **30x** reduction in this synthetic many-body query case.

Additional query samples on the same 4,096-body scene:

```text
batched sphere casts:     ~9.6 us/cast
small AABB overlap query: ~10.5 us/query
```

Pair-transition tracking is opt-in because it adds work proportional to the number of active body pairs. With `contactEvents=false`, the dense sphere/box baselines are effectively back at the v0.2 cost envelope despite the much larger v0.3 feature surface.

## Deliberate v0.3 limits

This is a compact game-oriented rigid-body kernel, not a claim of Rapier/Jolt feature parity.

Not implemented yet:

- convex hull / general convex support-map shapes;
- dynamic triangle meshes;
- slider, cone, six-DOF, gear, or ragdoll-specialized constraints;
- edge/edge multi-point manifolds;
- general dynamic-vs-dynamic TOI islands;
- moving-box/capsule/compound CCD and rotational TOI;
- character-controller policy;
- deterministic fixed-point simulation or cross-architecture bitwise determinism;
- solver parallelism across independent islands.

For M4J-style integration, the intended immediate shape is: fixed authoritative ticks, compound dynamic actors/ships, static mesh/heightfield world chunks, origin rebasing, bulk ECS transform exchange, batched spatial queries/sphere sweeps, and optional contact-pair event drains.
