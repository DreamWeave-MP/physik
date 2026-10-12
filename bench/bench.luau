--!strict
-- physik, timed from Luau: steps of piles, tumbles, compounds, terrain, chains and fast movers
-- through each broadphase, the spatial queries, the drains, and building shapes and worlds. Run
-- by `jess run --package . bench/bench.luau [filter] [results]` on a release build. Prints one
-- line per workload: its name, a tab, and the seconds one call takes, a call being what the name
-- says. Every workload's result is checked against what 0.3.2 computes, to the bit, so a path
-- that stopped doing the work, or does different work, fails here instead of posting a time.
-- With `results` as the second argument each workload's result is printed instead of timed,
-- which is how the table of expectations is made.
local physik = require '@dream/physik'

-- A build without the verified toolchain runs Luau slower than a release does, so a number from
-- one says nothing about what ships.
if not jess.bench.verified then
  error('physik bench: this jess was built without the verified toolchain: measure on a release build', 0)
end

local filter: string?, mode: string? = ...
local only = if filter and filter ~= '' then filter else nil
local printing = mode == 'results'

type World = typeof(physik.newWorld())
type Shape = typeof(physik.sphere(1))
type Random = (low: number, high: number) -> number

-- What every workload computes, as 0.3.2 computes it.
local Expected: { [string]: string } = {
  ['step: 1,000 spheres piled in the grid, 60 steps'] = '86747 contacts, 3541 candidates, 1601 manifolds, 1373 warm, 0 asleep, 54 islands, 0 swept, checksum -118484.44099515244',
  ['step: 1,000 spheres piled in the grid with contact events, 60 steps'] = '90964 events, 88217 contacts, 3549 candidates, 1625 manifolds, 1409 warm, 0 asleep, 54 islands, 0 swept, checksum -117784.94168266191',
  ['step: 250 boxes tumbling wider than a cell, 30 steps'] = '29301 contacts, 1549 candidates, 606 manifolds, 491 warm, 0 asleep, 1 islands, 0 swept, checksum -15663.445957111982',
  ['step: 250 boxes tumbling in the grid, 30 steps'] = '28089 contacts, 1497 candidates, 552 manifolds, 426 warm, 0 asleep, 1 islands, 0 swept, checksum -16647.894532551149',
  ['step: 600 spheres, capsules and boxes swept and pruned, 30 steps'] = '13181 contacts, 1286 candidates, 601 manifolds, 591 warm, 0 asleep, 142 islands, 0 swept, checksum -69309.090860039869',
  ['step: 2,000 spheres swept and pruned from scattered, 1 step'] = '150 contacts, 289 candidates, 150 manifolds, 0 warm, 0 asleep, 1853 islands, 0 swept, checksum 44388.679525341249',
  ['step: 256 compounds of three pieces, 30 steps'] = '5595 contacts, 1551 candidates, 170 manifolds, 153 warm, 0 asleep, 141 islands, 0 swept, checksum -36980.820001213622',
  ['step: 600 spheres, capsules and boxes on a heightfield 65 by 65, 30 steps'] = '28625 contacts, 473 candidates, 445 manifolds, 827 warm, 0 asleep, 567 islands, 0 swept, checksum -57145.373609047492',
  ['step: 600 spheres, capsules and boxes on a mesh of 8,192 triangles, 30 steps'] = '28570 contacts, 474 candidates, 438 manifolds, 812 warm, 0 asleep, 569 islands, 0 swept, checksum -58692.969448747368',
  ['step: 5,000 spheres drifting far apart with contact events, drained, 30 steps'] = '150000 drained, 107 contacts, 107 candidates, 40 manifolds, 25 warm, 0 asleep, 4960 islands, 0 swept, checksum 4163370.6571599795',
  ['step: 1,000 spheres asleep on the floor, 1 step'] = '1000 contacts, 1000 candidates, 1000 manifolds, 0 warm, 1000 asleep, 1000 islands, 0 swept, checksum 18256.543475389481',
  ['step: 60 chains of 8 links on distance, ball and hinge joints, 30 steps'] = '9000 contacts, 120 candidates, 120 manifolds, 0 warm, 0 asleep, 60 islands, 0 swept, checksum 33160.638126875238',
  ['step: 200 swept spheres fired at walls, spheres, compounds and terrain, 6 steps'] = '221 contacts, 63 candidates, 21 manifolds, 0 warm, 0 asleep, 199 islands, 32 swept, checksum 1846988.4500978959',
  ['query: 2,048 rays in a batch, as columns'] = '918 hits, checksum 67171.772320490592',
  ['query: 2,048 rays in a batch, as columns, unbounded and downward'] = '2048 hits, checksum 115321.69588919714',
  ['query: 2,048 rays in a batch, as records'] = '918 hits, checksum 28709.165104394498',
  ['query: 512 rays one by one'] = '222 hits, checksum 334669.87946955959',
  ['query: 1,024 sphere casts in a batch, as columns'] = '423 hits, checksum 13265.564990301822',
  ['query: 2,048 boxes of a few cells'] = '10145 found, checksum 20933567',
  ['query: 64 boxes of more cells than the grid is walked for'] = '177998 found, checksum 370614676',
  ['drain: transforms of 1,001 bodies'] = '-41144.848694401953',
  ['drain: linear state of 1,000 dynamic bodies, with user data and without handles'] = '1000 rows, checksum -717.93444389839874',
  ['drain: contacts of the pile, as columns'] = '1170 contacts, checksum 2102117.5334441387',
  ['drain: contacts of the pile, as records'] = '1170 contacts, checksum 491838.57829057064',
  ['drain: contact events of the pile, without the ones that stay'] = '1181 events, 50 that began or ended, checksum 55350',
  ['drain: stats of the pile'] = '1001 bodies, 0 asleep, 1170 contacts, 0 sensed',
  ['build: a mesh of 8,192 triangles from tables'] = '259.59986051982287 bounds, 97326863 chained, 0 overflowing',
  ['build: a mesh of 8,192 triangles from buffers'] = '259.59986051982287 bounds, 97326863 chained, 0 overflowing',
  ['build: a mesh of 1,408 triangles without indices'] = '259.59986051982287 bounds, 2311478 chained, 0 overflowing',
  ['build: a heightfield 65 by 65 from a table'] = '771.59986051982287 bounds, 0 chained, 0 overflowing',
  ['build: a heightfield 65 by 65 from a buffer'] = '771.59986051982287 bounds, 0 chained, 0 overflowing',
  ['build: add 5,000 bodies of every primitive, then remove every other one'] = '2500 left, checksum 31543.379452276655',
  ['sync: move the origin under 5,000 bodies, and back'] = '3934695.0420065979',
  ['sync: apply the transforms of 5,000 bodies'] = '5000 applied, checksum 3934695.0424996582',
}

local function median(samples: { number }): number
  table.sort(samples)

  return samples[(#samples + 1) // 2]
end

local function pin(value: number): string
  return if value % 1 == 0 then string.format('%d', value) else string.format('%.17g', value)
end

local function wanted(name: string): boolean
  return not only or string.find(name, only, 1, true) ~= nil
end

local function check(name: string, result: string)
  if result == Expected[name] then return end

  error(`physik bench: {name}: computed {result}, not {Expected[name]}`, 0)
end

--- Seconds per call of `body`: warmed up, then the median of seven batches of a size that takes
--- about a tenth of a second. The result of the first call and of one more after the last batch
--- is checked, outside the time, so a workload must compute the same thing every time it runs.
local function workload(name: string, body: () -> string)
  if not wanted(name) then return end

  if printing then
    print(`  ['{name}'] = '{body()}',`)
    return
  end

  jess.bench.collect()
  check(name, body())
  local start = os.clock()
  body()
  local once = math.max(os.clock() - start, 1e-7)
  local iterations = math.max(1, math.floor(0.1 / once))
  local samples = {}

  for batch = 1, 7 do
    jess.bench.collect()
    start = os.clock()

    for _ = 1, iterations do
      body()
    end

    samples[batch] = (os.clock() - start) / iterations
  end

  check(name, body())
  print(`{name}\t{string.format('%.9f', median(samples))}`)
end

--- A source of numbers in [low, high), the same ones for the same seed everywhere.
local function generator(seed: number): Random
  local state = seed

  return function(low: number, high: number): number
    state = bit32.band(state * 1664525 + 1013904223, 0xffffffff)

    return low + (high - low) * (state / 4294967296)
  end
end

--- A rotation with no trigonometry in it, so the same bits come out of every libm.
local function rotation(random: Random)
  return physik.quat(random(-1, 1), random(-1, 1), random(-1, 1), random(0.2, 1))
end

local function weigh(value: vector): number
  return value.X + 2 * value.Y + 3 * value.Z
end

local transforms = {}

--- Every live body's position, rotation and velocities, folded into one number.
local function checksum(world: World): number
  local state, count = world:writeTransforms(transforms)
  local positions, rotations, turns = state.positions, state.rotations, state.rotationW
  local velocities, spins = state.velocities, state.angularVelocities

  return sum[
    for index in range(1, count) =>
      (1 + index % 7) * (weigh(positions[index]) + weigh(rotations[index]) * 5 + turns[index] * 7)
        + (1 + index % 5) * (weigh(velocities[index]) * 11 + weigh(spins[index]) * 13)
  ]
end

--- What a world came to: the contacts its steps returned, what its counters say, and its state.
local function outcome(world: World, contacts: number): string
  local stats = world:stats()

  return `{contacts} contacts, {stats.candidates} candidates, {stats.manifolds} manifolds, {stats.warmStarted} warm, `
    .. `{stats.sleeping} asleep, {stats.islands} islands, {stats.ccdHits} swept, checksum {pin(checksum(world))}`
end

local function run(world: World, steps: number): string
  local contacts = 0

  for _ = 1, steps do
    contacts += world:step(1 / 60)
  end

  return outcome(world, contacts)
end

local function ground(world: World)
  world:addBody({
    type = 'static',
    shape = physik.box(vector.create(80, 0.5, 80)),
    position = vector.create(0, -0.5, 0),
    userData = 4096,
  })
end

--- `count` bodies in a block `side` wide and deep, `spacing` apart, resting on the floor or in it.
local function block(world: World, count: number, side: number, spacing: number, random: Random, shape: (index: number) -> Shape)
  for index = 0, count - 1 do
    local column, row, layer = index % side, index // side % side, index // (side * side)

    world:addBody({
      shape = shape(index),
      position = vector.create(
        (column - side / 2) * spacing + random(-0.05, 0.05),
        0.4 + layer * spacing,
        (row - side / 2) * spacing + random(-0.05, 0.05)
      ),
      rotation = rotation(random),
      angularVelocity = vector.create(random(-2, 2), random(-2, 2), random(-2, 2)),
      mass = random(0.5, 2),
    })
  end
end

local function mixed(random: Random): (index: number) -> Shape
  return function(index: number): Shape
    local kind = index % 3
    if kind == 0 then return physik.sphere(random(0.25, 0.4)) end
    if kind == 1 then return physik.capsule(random(0.15, 0.25), random(0.1, 0.3)) end

    return physik.box(vector.create(random(0.2, 0.35), random(0.2, 0.35), random(0.2, 0.35)))
  end
end

workload('step: 1,000 spheres piled in the grid, 60 steps', function(): string
  local random = generator(11)
  local world = physik.newWorld()
  ground(world)

  block(world, 1000, 10, 0.8, random, function(): Shape
    return physik.sphere(random(0.3, 0.45))
  end)

  return run(world, 60)
end)

workload('step: 1,000 spheres piled in the grid with contact events, 60 steps', function(): string
  local random = generator(12)
  local world = physik.newWorld({ contactEvents = true })
  ground(world)

  block(world, 1000, 10, 0.8, random, function(): Shape
    return physik.sphere(random(0.3, 0.45))
  end)

  local contacts, events, drain = 0, 0, {}

  for _ = 1, 60 do
    contacts += world:step(1 / 60)
    events += select(2, world:writeContactEvents(drain))
  end

  return `{events} events, {outcome(world, contacts)}`
end)

workload('step: 250 boxes tumbling wider than a cell, 30 steps', function(): string
  local random = generator(21)
  local world = physik.newWorld()
  ground(world)

  block(world, 250, 6, 1, random, function(): Shape
    return physik.box(vector.create(random(0.3, 0.55), random(0.3, 0.55), random(0.3, 0.55)))
  end)

  return run(world, 30)
end)

workload('step: 250 boxes tumbling in the grid, 30 steps', function(): string
  local random = generator(22)
  local world = physik.newWorld({ broadphaseCellSize = 2.5 })
  ground(world)

  block(world, 250, 6, 1, random, function(): Shape
    return physik.box(vector.create(random(0.3, 0.55), random(0.3, 0.55), random(0.3, 0.55)))
  end)

  return run(world, 30)
end)

workload('step: 600 spheres, capsules and boxes swept and pruned, 30 steps', function(): string
  local random = generator(31)
  local world = physik.newWorld({ broadphase = 'sap' })
  ground(world)
  block(world, 600, 9, 0.8, random, mixed(random))

  return run(world, 30)
end)

workload('step: 2,000 spheres swept and pruned from scattered, 1 step', function(): string
  local random = generator(32)
  local world = physik.newWorld({ broadphase = 'sap', gravity = vector.create(0, 0, 0) })

  for _ = 1, 2000 do
    world:addBody({
      shape = physik.sphere(random(0.2, 0.4)),
      position = vector.create(random(-400, 400), random(-2, 2), random(-2, 2)),
      mass = 1,
    })
  end

  return run(world, 1)
end)

workload('step: 256 compounds of three pieces, 30 steps', function(): string
  local random = generator(41)
  local world = physik.newWorld({ broadphaseCellSize = 3 })
  ground(world)

  block(world, 256, 7, 1.2, random, function(): Shape
    return physik.compound({
      { shape = physik.box(vector.create(0.5, 0.15, 0.25)) },
      { shape = physik.sphere(random(0.2, 0.3)), position = vector.create(0.55, 0.1, 0) },
      {
        shape = physik.capsule(0.12, random(0.15, 0.3)),
        position = vector.create(-0.5, 0.1, 0),
        rotation = rotation(random),
      },
    })
  end)

  return run(world, 30)
end)

-- One terrain, as a heightfield and as the triangle mesh of the same surface.
local TerrainSide = 65
local terrainHeights: { number } = {}

do
  local random = generator(51)

  for index = 1, TerrainSide * TerrainSide do
    terrainHeights[index] = random(0, 0.6)
  end
end

local terrainVertices: { vector } = [
  for index in range(0, TerrainSide * TerrainSide - 1) =>
    vector.create(index % TerrainSide - 32, terrainHeights[index + 1], index // TerrainSide - 32)
]

-- Two triangles to a cell, each corner counted from the cell's first.
local terrainIndices: { number } = [
  for row in range(0, TerrainSide - 2)
  for column in range(0, TerrainSide - 2)
  for corner in { 1, 2, TerrainSide + 2, 1, TerrainSide + 2, TerrainSide + 1 } =>
    row * TerrainSide + column + corner
]

local heightfield = physik.heightfield(TerrainSide, TerrainSide, terrainHeights, { origin = vector.create(-32, 0, -32) })
local mesh = physik.triangleMesh(terrainVertices, terrainIndices)

local function overTerrain(terrain: Shape, seed: number): World
  local random = generator(seed)
  local world = physik.newWorld({ broadphaseCellSize = 1.5 })
  world:addBody({ type = 'static', shape = terrain })
  block(world, 600, 24, 1.3, random, mixed(random))

  return world
end

workload('step: 600 spheres, capsules and boxes on a heightfield 65 by 65, 30 steps', function(): string
  return run(overTerrain(heightfield, 52), 30)
end)

workload('step: 600 spheres, capsules and boxes on a mesh of 8,192 triangles, 30 steps', function(): string
  return run(overTerrain(mesh, 52), 30)
end)

workload('step: 5,000 spheres drifting far apart with contact events, drained, 30 steps', function(): string
  local random = generator(61)
  local world = physik.newWorld({ gravity = vector.create(0, 0, 0), contactEvents = true, broadphaseCellSize = 2 })

  for index = 0, 4999 do
    world:addBody({
      shape = physik.sphere(0.5),
      position = vector.create(index % 17 * 4 + 1, index // 17 % 17 * 4 + 1, index // 289 * 4 + 1),
      velocity = vector.create(random(-4, 4), random(-4, 4), random(-4, 4)),
      mass = 1,
      userData = index,
    })
  end

  local contacts, drained, drain = 0, 0, {}

  for _ = 1, 30 do
    contacts += world:step(1 / 60)
    drained += select(2, world:writeTransforms(drain, { dynamicOnly = true }))
  end

  return `{drained} drained, {outcome(world, contacts)}`
end)

do
  local world = physik.newWorld()
  ground(world)

  for index = 0, 999 do
    world:addBody({
      shape = physik.sphere(0.4),
      position = vector.create(index % 32 - 16, 0.4, index // 32 - 16),
      mass = 1,
    })
  end

  for _ = 1, 90 do
    world:step(1 / 60)
  end

  workload('step: 1,000 spheres asleep on the floor, 1 step', function(): string
    return run(world, 1)
  end)
end

workload('step: 60 chains of 8 links on distance, ball and hinge joints, 30 steps', function(): string
  local random = generator(71)
  local world = physik.newWorld({ broadphaseCellSize = 2 })
  ground(world)

  for chain = 0, 59 do
    local column, row = chain % 8 * 3 - 12, chain // 8 * 3 - 12

    local above = world:addBody({
      type = 'static',
      shape = physik.sphere(0.2),
      position = vector.create(column, 4, row),
    })

    for link = 1, 8 do
      local body = world:addBody({
        shape = if link % 2 == 0 then physik.sphere(0.25) else physik.box(vector.create(0.2, 0.3, 0.2)),
        position = vector.create(column + link * 0.5, 4 - link * 0.6, row + random(-0.1, 0.1)),
        mass = 1,
      })

      local anchors = { localAnchorA = vector.create(0, -0.3, 0), localAnchorB = vector.create(0, 0.3, 0) }

      if link % 3 == 0 then
        world:addDistanceJoint(above, body, { length = 0.7 })
      elseif link % 3 == 1 then
        world:addBallJoint(above, body, anchors)
      else
        world:addHingeJoint(above, body, {
          localAnchorA = anchors.localAnchorA,
          localAnchorB = anchors.localAnchorB,
          localAxisA = vector.create(0, 0, 1),
          localAxisB = vector.create(0, 0, 1),
          motorSpeed = if chain % 2 == 0 then 1.5 else nil,
          maxMotorTorque = 4,
        })
      end

      above = body
    end
  end

  return run(world, 30)
end)

workload('step: 200 swept spheres fired at walls, spheres, compounds and terrain, 6 steps', function(): string
  local random = generator(81)
  local world = physik.newWorld({ broadphaseCellSize = 2 })
  world:addBody({ type = 'static', shape = heightfield, position = vector.create(0, -6, 0) })
  world:addBody({ type = 'static', shape = mesh, position = vector.create(0, 14, 0) })

  for index = 0, 59 do
    local position = vector.create(12 + index % 3, index // 3 % 5 * 2 - 4, index // 15 * 6 - 9)
    local kind = index % 3

    if kind == 0 then
      world:addBody({ type = 'static', shape = physik.box(vector.create(0.05, 1, 3)), position = position })
    elseif kind == 1 then
      world:addBody({ type = 'kinematic', shape = physik.sphere(0.8), position = position })
    else
      world:addBody({
        type = 'static',
        position = position,
        rotation = rotation(random),
        shape = physik.compound({
          { shape = physik.box(vector.create(0.4, 0.4, 0.4)) },
          { shape = physik.capsule(0.2, 0.5), position = vector.create(0, 0.9, 0) },
          { shape = physik.sphere(0.3), position = vector.create(0, -0.7, 0) },
        }),
      })
    end
  end

  for _ = 1, 200 do
    world:addBody({
      shape = physik.sphere(0.1),
      position = vector.create(random(-10, -2), random(-4, 4), random(-9, 9)),
      velocity = vector.create(random(200, 900), random(-400, 400), random(-60, 60)),
      mass = 0.2,
      ccd = true,
    })
  end

  return run(world, 6)
end)

-- A crowd for the queries: 4,096 small bodies of every primitive in a block, over one floor.
local crowd = physik.newWorld({ gravity = vector.create(0, 0, 0), broadphaseCellSize = 2 })

do
  local random = generator(91)
  local shape = mixed(random)
  ground(crowd)

  for index = 0, 4095 do
    crowd:addBody({
      shape = shape(index),
      position = vector.create(index % 16 * 2 - 16 + random(-0.3, 0.3), 1 + index // 256 * 2, index // 16 % 16 * 2 - 16 + random(-0.3, 0.3)),
      rotation = rotation(random),
      mass = 1,
      userData = index,
    })
  end

  crowd:step(1 / 60)
end

local QueryCount = 2048
local origins: { vector }, directions: { vector }, reaches: { number } = {}, {}, {}
local downward: { vector }, casts: { vector }, corners: { vector }, opposites: { vector } = {}, {}, {}, {}

do
  local random = generator(92)

  for index = 1, QueryCount do
    local origin = vector.create(random(-18, 18), random(0.5, 33), random(-18, 18))
    local direction = vector.create(random(-1, 1), random(-1, 1), random(-1, 1))
    local reach = vector.create(random(1, 2), random(1, 2), random(1, 2))

    origins[index], directions[index], reaches[index] = origin, direction, random(4, 40)
    downward[index] = vector.create(random(-0.2, 0.2), -1, random(-0.2, 0.2))
    casts[index] = origin + direction * random(2, 8)
    corners[index], opposites[index] = origin - reach, origin + reach
  end
end

local rays = {}

--- The hits of a batch of rays, drained as columns: how many, and what they hit and where.
local function struck(count: number): string
  local hits = #[for index in range(1, count) if rays.hit[index] => index]

  local total = sum[
    for index in range(1, count) if rays.hit[index] =>
      rays.distances[index] + weigh(rays.points[index]) + weigh(rays.normals[index]) * 3 + rays.userData[index] % 97
  ]

  return `{hits} hits, checksum {pin(total)}`
end

workload('query: 2,048 rays in a batch, as columns', function(): string
  return struck(select(2, crowd:raycastBatchSoA(origins, directions, reaches, nil, rays)))
end)

workload('query: 2,048 rays in a batch, as columns, unbounded and downward', function(): string
  return struck(select(2, crowd:raycastBatchSoA(origins, downward, nil, nil, rays)))
end)

local records = {}

workload('query: 2,048 rays in a batch, as records', function(): string
  local batch = crowd:raycastBatch(origins, directions, reaches, nil, records)
  local hits = #[for index in range(1, batch.count) if batch[index] => index]
  local total = sum[for index in range(1, batch.count) if batch[index] => batch[index].distance + weigh(batch[index].point)]

  return `{hits} hits, checksum {pin(total)}`
end)

workload('query: 512 rays one by one', function(): string
  local hits, total = 0, 0

  for index = 1, 512 do
    local hit = crowd:raycast(origins[index], directions[index], reaches[index])
    if not hit then continue end

    hits += 1
    total += hit.distance + weigh(hit.normal) + hit.handle % 1048576
  end

  return `{hits} hits, checksum {pin(total)}`
end)

local sweeps = {}

workload('query: 1,024 sphere casts in a batch, as columns', function(): string
  local batch, count = crowd:sphereCastBatchSoA(origins[1:1024], casts[1:1024], 0.25, nil, nil, sweeps)
  local hits = #[for index in range(1, count) if batch.hit[index] => index]

  local total = sum[
    for index in range(1, count) if batch.hit[index] =>
      batch.fractions[index] + batch.distances[index] + weigh(batch.centers[index]) + weigh(batch.normals[index]) * 3
  ]

  return `{hits} hits, checksum {pin(total)}`
end)

local overlaps = {}

workload('query: 2,048 boxes of a few cells', function(): string
  local found, total = 0, 0

  for index = 1, QueryCount do
    local handles, count = crowd:queryAabb(corners[index], opposites[index], nil, overlaps)

    found += count
    total += sum[for handle in handles => handle % 1048576]
  end

  return `{found} found, checksum {total}`
end)

workload('query: 64 boxes of more cells than the grid is walked for', function(): string
  local found, total = 0, 0

  for index = 1, 64 do
    local reach = vector.create(22, 22, 22)
    local handles, count = crowd:queryAabb(origins[index] - reach, origins[index] + reach, nil, overlaps)

    found += count
    total += sum[for handle in handles => handle % 1048576]
  end

  return `{found} found, checksum {total}`
end)

-- A pile mid-fall for the drains: every body awake, some three thousand contacts.
local pile = physik.newWorld({ contactEvents = true })

do
  local random = generator(101)
  ground(pile)

  block(pile, 1000, 10, 0.8, random, function(): Shape
    return physik.sphere(random(0.3, 0.45))
  end)

  for _ = 1, 10 do
    pile:step(1 / 60)
  end
end

local drained = {}

workload('drain: transforms of 1,001 bodies', function(): string
  return pin(checksum(pile))
end)

workload('drain: linear state of 1,000 dynamic bodies, with user data and without handles', function(): string
  local state, count = pile:writeLinearState(drained, { dynamicOnly = true, userData = true, handles = false })
  local total = sum[for index in range(1, count) => weigh(state.positions[index]) + weigh(state.velocities[index]) * 3]

  return `{count} rows, checksum {pin(total)}`
end)

local touching = {}

workload('drain: contacts of the pile, as columns', function(): string
  local state, count = pile:writeContacts(touching)

  local total = sum[
    for index in range(1, count) =>
      state.a[index] % 1048576 + state.b[index] % 1048576 * 3 + weigh(state.normals[index]) + weigh(state.points[index]) + state.penetration[index] + state.normalImpulse[index]
  ]

  return `{count} contacts, checksum {pin(total)}`
end)

workload('drain: contacts of the pile, as records', function(): string
  local contacts = pile:contacts()
  local total = sum[for contact in contacts => contact.a % 1048576 + weigh(contact.point) + contact.penetration + contact.normalImpulse]

  return `{#contacts} contacts, checksum {pin(total)}`
end)

local happened = {}

workload('drain: contact events of the pile, without the ones that stay', function(): string
  local every = select(2, pile:writeContactEvents(happened))
  local state, count = pile:writeContactEvents(happened, false)
  local total = sum[for index in range(1, count) => state.a[index] % 1048576 + state.b[index] % 1048576 * 3 + state.kind[index]]

  return `{every} events, {count} that began or ended, checksum {total}`
end)

workload('drain: stats of the pile', function(): string
  local stats = pile:stats()

  return `{stats.bodies} bodies, {stats.sleeping} asleep, {stats.contacts} contacts, {stats.sensorContacts} sensed`
end)

--- What a shape's accelerator holds, by walking it the way a query does.
local function measured(shape): string
  local low, high = shape.localMin, shape.localMax
  local chained = 0

  for _, head in shape.heads or {} do
    local node = head

    while node ~= 0 do
      chained += shape.nodeTri[node] + 1
      node = shape.nextNode[node]
    end
  end

  return `{pin(weigh(low) + weigh(high) * 3)} bounds, {chained} chained, {#(shape.overflow or {})} overflowing`
end

workload('build: a mesh of 8,192 triangles from tables', function(): string
  return measured(physik.triangleMesh(terrainVertices, terrainIndices))
end)

workload('build: a mesh of 8,192 triangles from buffers', function(): string
  return measured(physik.triangleMesh(mesh.vertices, mesh.indices))
end)

workload('build: a mesh of 1,408 triangles without indices', function(): string
  return measured(physik.triangleMesh(terrainVertices[1:4224]))
end)

workload('build: a heightfield 65 by 65 from a table', function(): string
  return measured(physik.heightfield(TerrainSide, TerrainSide, terrainHeights))
end)

workload('build: a heightfield 65 by 65 from a buffer', function(): string
  return measured(physik.heightfield(TerrainSide, TerrainSide, heightfield.heights))
end)

workload('build: add 5,000 bodies of every primitive, then remove every other one', function(): string
  local random = generator(111)
  local shape = mixed(random)
  local world = physik.newWorld()

  local handles = [
    for index in range(0, 4999) =>
      world:addBody({
        shape = shape(index),
        position = vector.create(random(-50, 50), random(-50, 50), random(-50, 50)),
        rotation = rotation(random),
      })
  ]

  for index = 1, 5000, 2 do
    world:removeBody(handles[index])
  end

  return `{world:stats().bodies} left, checksum {pin(checksum(world))}`
end)

do
  local random = generator(121)
  local world = physik.newWorld({ gravity = vector.create(0, 0, 0), broadphaseCellSize = 2 })
  local shape = mixed(random)

  for index = 0, 4999 do
    world:addBody({
      shape = shape(index),
      position = vector.create(index % 17 * 4, index // 17 % 17 * 4, index // 289 * 4),
      rotation = rotation(random),
      userData = index,
    })
  end

  -- Whole coordinates, so a shift and its opposite leave every position as it was, to the bit.
  workload('sync: move the origin under 5,000 bodies, and back', function(): string
    world:shiftOrigin(vector.create(64, -32, 16))
    world:shiftOrigin(vector.create(-64, 32, -16))

    return pin(checksum(world))
  end)

  -- Applying normalizes a rotation again, so every call applies the one state written here.
  local written = world:writeTransforms({})

  workload('sync: apply the transforms of 5,000 bodies', function(): string
    return `{world:applyTransforms(written, false)} applied, checksum {pin(checksum(world))}`
  end)
end
