local physik = require('@dream/physik')
local test = require('@jess/test')
local expect = test.expect
local suite = test.new()

local Generation = 1048576

--- Fails unless `actual` is within `tolerance` of `expected`: positions and distances come out
--- of 32-bit vector arithmetic, and a case says how far that may stray.
local function near(actual: number, expected: number, tolerance: number?)
  local allowed = tolerance or 1e-4
  if math.abs(actual - expected) <= allowed then return end

  error(`expected {actual} to be within {allowed} of {expected}`, 2)
end

local function nearVector(actual: vector, expected: vector, tolerance: number?)
  near(actual.X, expected.X, tolerance)
  near(actual.Y, expected.Y, tolerance)
  near(actual.Z, expected.Z, tolerance)
end

local function ground(world)
  return world:addBody({
    type = 'static',
    shape = physik.box(vector.create(20, 0.5, 20)),
    position = vector.create(0, -0.5, 0),
    userData = 'ground',
  })
end

--- The user data of the bodies behind `count` handles, in order.
local function named(world, handles: { number }, count: number): { string }
  local names = [for handle in handles[1:count] => world:getBody(handle).userData]
  table.sort(names)

  return names
end

suite.describe('shapes', function()
  suite.it('refuses sizes that are not positive', function()
    expect(function()
      physik.sphere(0)
    end).toThrow('physik.sphere: radius must be > 0')

    expect(function()
      physik.box(vector.create(1, 0, 1))
    end).toThrow('physik.box: extents must be > 0')

    expect(function()
      physik.capsule(1, -1)
    end).toThrow('physik.capsule: halfHeight must be >= 0')
  end)

  suite.it('builds the same mesh from tables and from packed buffers', function()
    local vertices = { vector.create(0, 0, 0), vector.create(4, 0, 0), vector.create(4, 1, 4), vector.create(0, 2, 4) }
    local fromTables = physik.triangleMesh(vertices, { 1, 2, 3, 1, 3, 4 })
    local fromBuffers = physik.triangleMesh(fromTables.vertices, fromTables.indices)

    for _, mesh in { fromTables, fromBuffers } do
      expect(mesh.vertexCount).toBe(4)
      expect(mesh.triangleCount).toBe(2)
      expect(mesh.localMin).toBe(vector.create(0, 0, 0))
      expect(mesh.localMax).toBe(vector.create(4, 2, 4))
      expect(mesh.cellSize).toBe(4)
      expect(mesh.overflow).toEqual({})
    end

    expect(buffer.readu32(fromTables.indices, 8)).toBe(2)
    expect(buffer.readf32(fromTables.vertices, 28)).toBe(1)
    expect(fromBuffers.heads).toEqual(fromTables.heads)
    expect(fromBuffers.nextNode).toEqual(fromTables.nextNode)
    expect(fromBuffers.nodeTri).toEqual(fromTables.nodeTri)
    expect(#fromTables.nodeTri).toBe(8)
  end)

  suite.it('numbers the triangles of a mesh without indices in vertex order', function()
    local mesh = physik.triangleMesh({
      vector.create(0, 0, 0),
      vector.create(1, 0, 0),
      vector.create(0, 0, 1),
      vector.create(5, 3, 5),
      vector.create(6, 3, 5),
      vector.create(5, 3, 6),
    })

    expect(mesh.triangleCount).toBe(2)
    expect(buffer.len(mesh.indices)).toBe(24)
    expect([for index in mesh.indices[1:6, 'u32'] => index]).toEqual({ 0, 1, 2, 3, 4, 5 })
    expect(mesh.localMin).toBe(vector.create(0, 0, 0))
    expect(mesh.localMax).toBe(vector.create(6, 3, 6))
  end)

  suite.it('keeps a triangle that covers too many cells in the overflow list', function()
    local vertices = { vector.create(0, 0, 0), vector.create(8, 0, 0), vector.create(0, 0, 8), vector.create(0.1, 0, 0.1) }
    local mesh = physik.triangleMesh(vertices, { 1, 2, 3, 1, 4, 2 }, { cellSize = 1, maxCellsPerTriangle = 16 })

    expect(mesh.overflow).toEqual({ 0 })
    expect(#mesh.nodeTri).toBe(9)
  end)

  suite.it('leaves a vertex that is not a number out of a mesh\'s bounds', function()
    local mesh = physik.triangleMesh({
      vector.create(0 / 0, 1, 1),
      vector.create(2, 0 / 0, 3),
      vector.create(-1, 5, 0 / 0),
    })

    expect(mesh.localMin).toBe(vector.create(-1, 1, 1))
    expect(mesh.localMax).toBe(vector.create(2, 5, 3))
  end)

  suite.it('refuses meshes it cannot index', function()
    local triangle = { vector.create(0, 0, 0), vector.create(1, 0, 0), vector.create(0, 0, 1) }

    expect(function()
      physik.triangleMesh({ triangle[1], triangle[2], triangle[3], triangle[1] })
    end).toThrow('physik.triangleMesh: without indices, vertex count must be divisible by 3')

    expect(function()
      physik.triangleMesh(triangle, { 1, 2, 4 })
    end).toThrow('physik.triangleMesh: table indices are 1-based and must reference a vertex')

    expect(function()
      physik.triangleMesh(triangle, { 0, 1, 2 })
    end).toThrow('physik.triangleMesh: table indices are 1-based and must reference a vertex')

    expect(function()
      physik.triangleMesh({ triangle[1], triangle[2], 3 })
    end).toThrow('physik.triangleMesh: every vertex must be a vector')

    expect(function()
      physik.triangleMesh(triangle, {})
    end).toThrow('physik.triangleMesh: mesh must contain at least one triangle')

    expect(function()
      physik.triangleMesh(triangle, buffer.create(8))
    end).toThrow('physik.triangleMesh: index buffer length must be divisible by 3 u32 indices')

    expect(function()
      physik.triangleMesh(buffer.create(10))
    end).toThrow('physik.triangleMesh: vertex buffer must contain packed f32 xyz triples')
  end)

  suite.it('bounds a heightfield by its lowest and highest samples', function()
    local options = { cellSize = 2, origin = vector.create(10, 100, -10) }
    local fromTable = physik.heightfield(3, 2, { 1, 5, -2, 0.5, 3, 2 }, options)
    local fromBuffer = physik.heightfield(3, 2, fromTable.heights, options)

    for _, field in { fromTable, fromBuffer } do
      expect(field.localMin).toBe(vector.create(10, 98, -10))
      expect(field.localMax).toBe(vector.create(14, 105, -8))
      expect(field.cellSizeX).toBe(2)
      expect(field.cellSizeZ).toBe(2)
    end

    expect([for height in fromTable.heights[1:6, 'f32'] => height]).toEqual({ 1, 5, -2, 0.5, 3, 2 })
    expect(fromBuffer.heights).toBe(fromTable.heights)
  end)

  suite.it('leaves a sample that is not a number out of a heightfield\'s bounds', function()
    local first = physik.heightfield(2, 2, { 0 / 0, 1, 2, 3 })
    local last = physik.heightfield(2, 2, { 4, -1, 2, 0 / 0 })
    local none = physik.heightfield(2, 2, { 0 / 0, 0 / 0, 0 / 0, 0 / 0 })
    local unbounded = physik.heightfield(2, 2, { math.huge, -math.huge, 0, 1 })

    expect(first.localMin.Y).toBe(1)
    expect(first.localMax.Y).toBe(3)
    expect(last.localMin.Y).toBe(-1)
    expect(last.localMax.Y).toBe(4)
    expect(none.localMin.Y).toBe(math.huge)
    expect(none.localMax.Y).toBe(-math.huge)
    expect(unbounded.localMin.Y).toBe(-math.huge)
    expect(unbounded.localMax.Y).toBe(math.huge)
  end)

  suite.it('refuses heightfields that are not a grid of numbers', function()
    expect(function()
      physik.heightfield(1, 2, { 1, 2 })
    end).toThrow('physik.heightfield: width must be an integer >= 2')

    expect(function()
      physik.heightfield(2, 2, { 1, 2, 3 })
    end).toThrow('physik.heightfield: heights must contain width*depth values')

    expect(function()
      physik.heightfield(2, 2, { 1, 2, 'three', 4 })
    end).toThrow('physik.heightfield: heights must be numbers')

    expect(function()
      physik.heightfield(2, 2, buffer.create(12))
    end).toThrow('physik.heightfield: height buffer must contain width*depth f32 values')

    expect(function()
      physik.heightfield(2, 2, { 1, 2, 3, 4 }, { cellSize = 0 })
    end).toThrow('physik.heightfield: cell sizes must be > 0')
  end)

  suite.it('bounds a compound by its children where they are placed', function()
    local compound = physik.compound({
      { shape = physik.sphere(1), position = vector.create(3, 0, 0) },
      { shape = physik.box(vector.create(0.5, 2, 0.5)) },
      { shape = physik.capsule(0.25, 1), position = vector.create(0, 0, -4) },
    })

    expect(compound.count).toBe(3)
    expect(compound.localMin).toBe(vector.create(-0.5, -2, -4.25))
    expect(compound.localMax).toBe(vector.create(4, 2, 1))

    expect(function()
      physik.compound({})
    end).toThrow('physik.compound: children must be a non-empty array')

    expect(function()
      physik.compound({ { shape = compound } })
    end).toThrow('physik.compound: children must be sphere, box, or capsule shapes')
  end)
end)

suite.describe('bodies', function()
  suite.it('hands out handles that go stale when their slot is reused', function()
    local world = physik.newWorld()
    local first = world:addBody({ shape = physik.sphere(1), userData = 'first' })

    expect(first).toBe(Generation + 1)
    expect(world:getBody(first).userData).toBe('first')
    world:removeBody(first)

    expect(function()
      world:getBody(first)
    end).toThrow('physik: stale or invalid body handle')

    expect(world:addBody({})).toBe(2 * Generation + 1)
    expect(world:addBody({})).toBe(Generation + 2)
    expect(world:stats().bodies).toBe(2)
  end)

  suite.it('describes a body as it was added', function()
    local world = physik.newWorld()

    local handle = world:addBody({
      type = 'kinematic',
      shape = physik.box(vector.create(1, 2, 3)),
      position = vector.create(4, 5, 6),
      rotation = physik.quat(0, 0, 3, 4),
      velocity = vector.create(1, 0, 0),
      sensor = true,
    })

    local body = world:getBody(handle)

    expect(body.handle).toBe(handle)
    expect(body.type).toBe('kinematic')
    expect(body.position).toBe(vector.create(4, 5, 6))
    nearVector(body.rotation.v, vector.create(0, 0, 0.6), 1e-6)
    near(body.rotation.w, 0.8, 1e-6)
    expect(body.velocity).toBe(vector.create(1, 0, 0))
    expect(body.angularVelocity).toBe(vector.create(0, 0, 0))
    expect(body.sleeping).toBe(false)
    expect(body.sensor).toBe(true)
  end)

  suite.it('refuses bodies it cannot simulate', function()
    local world = physik.newWorld()
    local mesh = physik.triangleMesh({ vector.create(0, 0, 0), vector.create(1, 0, 0), vector.create(0, 0, 1) })

    expect(function()
      world:addBody({ type = 'floating' })
    end).toThrow("physik: body type must be 'dynamic', 'static', or 'kinematic'")

    expect(function()
      world:addBody({ shape = mesh })
    end).toThrow('physik.addBody: triangle meshes and heightfields must be static or kinematic')

    expect(function()
      world:addBody({ shape = {} })
    end).toThrow('physik.addBody: invalid shape')

    expect(function()
      world:addBody({ rotation = vector.create(0, 0, 1) })
    end).toThrow('physik.addBody: rotation must come from physik.quat/axisAngle')

    -- A body that is refused leaves nothing behind, and the world goes on stepping.
    expect(world:stats().bodies).toBe(0)
    expect(world:step(1 / 60)).toBe(0)
    expect(world:addBody({})).toBe(Generation + 1)
  end)
end)

suite.describe('stepping', function()
  suite.it('integrates a falling body one step at a time', function()
    local world = physik.newWorld({ gravity = vector.create(0, -10, 0) })
    local handle = world:addBody({ shape = physik.sphere(0.5), linearDamping = 0, angularDamping = 0 })

    for _ = 1, 60 do
      expect(world:step(1 / 60)).toBe(0)
    end

    local body = world:getBody(handle)

    near(body.velocity.Y, -10, 1e-3)
    near(body.position.Y, -10 * 1830 / 3600, 1e-3)
    expect(world:stats().frame).toBe(60)
  end)

  suite.it('takes substeps and refuses a step that is not forward', function()
    local world = physik.newWorld({ gravity = vector.create(0, -10, 0) })
    local handle = world:addBody({ shape = physik.sphere(0.5), linearDamping = 0 })
    world:step(1 / 30, 2)

    near(world:getBody(handle).position.Y, -10 * 3 / 3600, 1e-5)
    expect(world:stats().frame).toBe(1)

    expect(function()
      world:step(0)
    end).toThrow('physik.step: dt must be > 0')

    expect(function()
      world:step(1 / 60, 1.5)
    end).toThrow('physik.step: substeps must be a positive integer')
  end)

  suite.it('catches up in fixed steps and drops what it cannot', function()
    local world = physik.newWorld({ fixedDt = 0.01, maxCatchUpSteps = 4 })
    local steps, alpha, dropped = world:advance(0.025)

    expect(steps).toBe(2)
    near(alpha, 0.5, 1e-9)
    expect(dropped).toBe(false)
    steps, alpha, dropped = world:advance(0.1)
    expect(steps).toBe(4)
    expect(dropped).toBe(true)
    near(alpha, 0.5, 1e-6)
    expect(world:stats().frame).toBe(6)
  end)

  for _, broadphase in { 'grid', 'sap' } do
    suite.it(`rests a sphere, a box and a capsule on the floor and puts them to sleep, by {broadphase}`, function()
      local world = physik.newWorld({ broadphase = broadphase })
      ground(world)

      local sphere = world:addBody({ shape = physik.sphere(0.5), position = vector.create(-3, 0.6, 0) })
      local box = world:addBody({ shape = physik.box(vector.create(0.5, 0.25, 0.5)), position = vector.create(0, 0.3, 0) })

      local capsule = world:addBody({
        shape = physik.capsule(0.25, 0.5),
        position = vector.create(3, 0.9, 0),
      })

      for _ = 1, 240 do
        world:step(1 / 60)
      end

      local stats = world:stats()

      -- Each settles a little into the floor: the solver leaves the slop it is allowed.
      near(world:getBody(sphere).position.Y, 0.4773, 1e-3)
      near(world:getBody(box).position.Y, 0.2441, 1e-3)
      near(world:getBody(capsule).position.Y, 0.7273, 1e-3)
      near(world:getBody(sphere).position.X, -3, 1e-3)
      expect(stats.sleeping).toBe(3)
      expect(stats.activeContacts).toBe(0)
      expect(stats.bodies).toBe(4)
      expect(stats.candidates).toBe(3)
      expect(stats.manifolds).toBe(3)
      expect(stats.contacts).toBe(6)
      expect(world:getBody(box).sleeping).toBe(true)

      world:applyImpulse(box, vector.create(0, 3, 0))
      expect(world:getBody(box).sleeping).toBe(false)
      world:step(1 / 60)
      expect(world:stats().sleeping).toBe(2)
      near(world:getBody(box).position.Y, 0.2914, 1e-3)
    end)

    suite.it(`stacks boxes and keeps every pair of neighbours in contact, by {broadphase}`, function()
      local world = physik.newWorld({ broadphase = broadphase })
      ground(world)

      local boxes = [
        for level in range(1, 2) =>
          world:addBody({
            shape = physik.box(vector.create(0.5, 0.5, 0.5)),
            position = vector.create(0, level - 0.5, 0),
          })
      ]

      for _ = 1, 30 do
        world:step(1 / 60)
      end

      local contacts = world:contacts()
      local stats = world:stats()

      for level, box in boxes do
        local body = world:getBody(box)

        near(body.position.Y, level - 0.5, 0.08)
        near(body.position.X, 0, 0.02)
        near(body.position.Z, 0, 0.02)
      end

      expect(stats.manifolds).toBe(2)
      expect(stats.contacts).toBe(8)
      expect(#contacts).toBe(8)
      expect(stats.islands).toBe(1)

      for _, contact in contacts do
        nearVector(contact.normal, vector.create(0, 1, 0), 0.02)
        expect(contact.sensor).toBe(false)
      end
    end)

    suite.it(`pairs bodies added out of order and far apart, by {broadphase}`, function()
      local world = physik.newWorld({ broadphase = broadphase, gravity = vector.create(0, 0, 0) })
      local handles = {}

      for index = 40, 1, -1 do
        local column = (index * 7 % 40) * 3

        handles[#handles + 1] = world:addBody({ shape = physik.sphere(0.5), position = vector.create(column, 0, 0) })
        handles[#handles + 1] = world:addBody({ shape = physik.box(vector.create(0.4, 0.4, 0.4)), position = vector.create(column + 0.8, 0, 0) })
      end

      expect(world:step(1 / 60)).toBe(40)

      local pairs = [for contact in world:contacts() => math.abs(contact.a - contact.b)]
      expect(pairs).toEqual(table.create(40, 1))

      -- Enough removals that the sweep's order is compacted before it is sorted again.
      for index = 1, 70 do
        world:removeBody(handles[index])
      end

      world:step(1 / 60)
      expect(world:stats().bodies).toBe(10)
      expect(world:stats().contacts).toBe(5)
    end)

    suite.it(`keeps a sphere on the terrain it rolls across, by {broadphase}`, function()
      local flat = table.create(33 * 33, 0)
      local vertices = { vector.create(-16, 0, -16), vector.create(16, 0, -16), vector.create(16, 0, 16), vector.create(-16, 0, 16) }

      local terrains = {
        physik.heightfield(33, 33, flat, { origin = vector.create(-16, 0, -16) }),
        physik.triangleMesh(vertices, { 1, 3, 2, 1, 4, 3 }, { cellSize = 1, maxCellsPerTriangle = 4096 }),
      }

      for _, terrain in terrains do
        local world = physik.newWorld({ broadphase = broadphase })
        world:addBody({ type = 'static', shape = terrain })

        local ball = world:addBody({
          shape = physik.sphere(0.4),
          position = vector.create(-10, 0.4, 0.3),
          velocity = vector.create(6, 0, 0),
          linearDamping = 0,
          friction = 0,
        })

        local lowest = math.huge

        for _ = 1, 150 do
          world:step(1 / 60)
          lowest = math.min(lowest, world:getBody(ball).position.Y)
        end

        -- Ten units on, far from the cells it started over, it is still resting on the surface.
        expect(lowest > 0.3).toBe(true)
        expect(world:getBody(ball).position.X > 0).toBe(true)
        expect(world:stats().contacts > 0).toBe(true)
      end
    end)

    suite.it(`pairs bodies that share a least X in the order they were added, by {broadphase}`, function()
      local world = physik.newWorld({ broadphase = broadphase, gravity = vector.create(0, 0, 0), broadphaseCellSize = 2 })

      -- Six columns of twelve spheres, each touching the one above; a column's spheres all
      -- start at the same X, and the columns are added out of order, in two passes.
      for pass = 0, 1 do
        for index = 0, 35 do
          local column, level = index * 5 % 6, index // 6 * 2 + pass

          world:addBody({
            shape = physik.sphere(0.5),
            position = vector.create(column * 3, level * 0.9, 0),
          })
        end
      end

      expect(world:step(1 / 60)).toBe(66)

      local contacts = world:contacts()
      local order = sum[for index, contact in enumerate(contacts) => index * (contact.a % Generation * 100 + contact.b % Generation)]

      expect(order).toBe(if broadphase == 'sap' then 4678762 else 5070681)
    end)
  end

  suite.it('collides every pair of primitives', function()
    local shapes = {
      sphere = physik.sphere(0.5),
      box = physik.box(vector.create(0.5, 0.5, 0.5)),
      capsule = physik.capsule(0.3, 0.4),
      compound = physik.compound({
        { shape = physik.sphere(0.3), position = vector.create(-0.2, 0, 0) },
        { shape = physik.box(vector.create(0.3, 0.3, 0.3)), position = vector.create(0.2, 0, 0) },
      }),
    }

    local expected = {
      ['sphere sphere'] = 1,
      ['sphere box'] = 1,
      ['sphere capsule'] = 1,
      ['sphere compound'] = 2,
      ['box sphere'] = 1,
      ['box box'] = 4,
      ['box capsule'] = 1,
      ['box compound'] = 2,
      ['capsule sphere'] = 1,
      ['capsule box'] = 1,
      ['capsule capsule'] = 1,
      ['capsule compound'] = 1,
      ['compound sphere'] = 2,
      ['compound box'] = 2,
      ['compound capsule'] = 1,
      ['compound compound'] = 3,
    }

    for _, left in { 'sphere', 'box', 'capsule', 'compound' } do
      for _, right in { 'sphere', 'box', 'capsule', 'compound' } do
        local world = physik.newWorld({ gravity = vector.create(0, 0, 0), broadphaseCellSize = 4 })

        local first = world:addBody({ shape = shapes[left], position = vector.create(0, 0, 0) })
        local second = world:addBody({ shape = shapes[right], position = vector.create(0.5, 0.1, 0.05) })
        local count = world:step(1 / 60)
        local contact = world:contacts()[1]

        if count ~= expected[`{left} {right}`] then
          error(`{left} against {right} made {count} contacts, not {expected[`{left} {right}`]}`)
        end

        expect(contact.a).toBe(first)
        expect(contact.b).toBe(second)
        expect(contact.normal.X > 0.9).toBe(true)
        expect(contact.penetration > 0).toBe(true)
      end
    end
  end)

  suite.it('rests primitives on a heightfield and on the mesh of the same surface', function()
    local heights = table.create(25, 1)
    local vertices = [for index in range(0, 24) => vector.create(index % 5 * 2 - 4, 1, index // 5 * 2 - 4)]

    local indices = [
      for row in range(0, 3)
      for column in range(0, 3)
      for corner in { 1, 7, 2, 1, 6, 7 } =>
        row * 5 + column + corner
    ]

    local terrains = {
      physik.heightfield(5, 5, heights, { cellSize = 2, origin = vector.create(-4, 0, -4) }),
      physik.triangleMesh(vertices, indices),
    }

    for _, terrain in terrains do
      local world = physik.newWorld()
      world:addBody({ type = 'static', shape = terrain })

      local sphere = world:addBody({ shape = physik.sphere(0.5), position = vector.create(-2.5, 1.6, 0.3) })
      local box = world:addBody({ shape = physik.box(vector.create(0.4, 0.4, 0.4)), position = vector.create(0.2, 1.5, 0.3) })
      local capsule = world:addBody({ shape = physik.capsule(0.3, 0.3), position = vector.create(2.6, 1.7, -0.4) })

      for _ = 1, 180 do
        world:step(1 / 60)
      end

      near(world:getBody(sphere).position.Y, 1.4773, 1e-3)
      near(world:getBody(capsule).position.Y, 1.5773, 1e-3)
      nearVector(world:getBody(sphere).position, vector.create(-2.5, 1.4773, 0.3), 1e-3)

      -- A box does not come to rest on triangles: it rides a little above them, and wanders.
      near(world:getBody(box).position.Y, 1.6, 0.1)
      expect(world:stats().sleeping).toBe(2)
      expect(world:stats().contacts).toBe(3)
    end
  end)

  suite.it('stops a small sphere at the first of more oversized bodies than are paired one by one', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0), broadphaseCellSize = 1 })

    -- Twenty walls, each wider than a cell: more than the sixteen that are tested against
    -- every body directly, so the pairs come from the sweep along X.
    for index = 1, 20 do
      world:addBody({ type = 'static', shape = physik.box(vector.create(1, 1, 1)), position = vector.create(100 + index * 5, 0, 0) })
    end

    local ball = world:addBody({
      shape = physik.sphere(0.3),
      velocity = vector.create(20, 0, 0),
      linearDamping = 0,
    })

    local contacts = 0

    for _ = 1, 400 do
      contacts += world:step(1 / 60)
    end

    expect(contacts > 0).toBe(true)
    expect(world:getBody(ball).position.X < 104).toBe(true)
  end)

  suite.it('keeps a sensor out of the solve and reports what it touches', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0), contactEvents = true })
    local sensor = world:addBody({ type = 'static', shape = physik.sphere(1), sensor = true })

    local mover = world:addBody({
      shape = physik.sphere(0.5),
      position = vector.create(-3, 0, 0),
      velocity = vector.create(6, 0, 0),
      linearDamping = 0,
    })

    local kinds = {}

    for _ = 1, 60 do
      world:step(1 / 60)

      local events, count = world:writeContactEvents()
      if count == 0 then continue end

      expect(count).toBe(1)
      expect(events.a[1]).toBe(sensor)
      expect(events.b[1]).toBe(mover)
      expect(events.sensor[1]).toBe(true)
      table.insert(kinds, events.kind[1])
    end

    local event = physik.contactEvent

    expect(kinds[1]).toBe(event.enter)
    expect(kinds[#kinds]).toBe(event.exit)
    expect(#[for kind in kinds if kind == event.stay => kind]).toBe(#kinds - 2)
    expect(#kinds).toBe(32)
    near(world:getBody(mover).velocity.X, 6, 1e-4)
    near(world:getBody(mover).position.X, 3, 1e-3)
    expect(select(2, world:writeContactEvents({}, false))).toBe(0)
  end)

  suite.it('filters pairs by group and mask', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0) })
    world:addBody({ shape = physik.sphere(0.5), group = 2, mask = 4 })
    world:addBody({ shape = physik.sphere(0.5), position = vector.create(0.5, 0, 0), group = 4, mask = 1 })
    expect(world:step(1 / 60)).toBe(0)

    world:addBody({ shape = physik.sphere(0.5), position = vector.create(0, 0.5, 0), group = 4, mask = 2 })
    expect(world:step(1 / 60)).toBe(1)
  end)

  suite.it('stops a swept sphere at a wall a discrete step would cross', function()
    local function fire(swept: boolean): number
      local world = physik.newWorld({ gravity = vector.create(0, 0, 0) })
      world:addBody({ type = 'static', shape = physik.box(vector.create(0.05, 2, 2)), position = vector.create(5, 0, 0) })

      local bullet = world:addBody({
        shape = physik.sphere(0.1),
        velocity = vector.create(1000, 0, 0),
        linearDamping = 0,
        ccd = swept,
      })

      world:step(1 / 60)
      expect(world:stats().ccdHits).toBe(if swept then 1 else 0)

      return world:getBody(bullet).position.X
    end

    near(fire(false), 1000 / 60, 1e-3)
    near(fire(true), 4.85, 0.01)
  end)

  suite.it('holds bodies together with joints until they are removed', function()
    local world = physik.newWorld({ gravity = vector.create(0, -10, 0) })
    local anchor = world:addBody({ type = 'static', shape = physik.sphere(0.2), position = vector.create(0, 5, 0) })
    local weight = world:addBody({ shape = physik.sphere(0.3), position = vector.create(0, 4, 0) })
    local joint = world:addDistanceJoint(anchor, weight, { length = 1, stiffness = 0.5 })

    expect(joint).toBe(1)
    expect(world:stats().joints).toBe(1)

    for _ = 1, 120 do
      world:step(1 / 60)
    end

    near(world:getBody(weight).position.Y, 4, 0.1)
    world:removeJoint(joint)
    expect(world:stats().joints).toBe(0)

    for _ = 1, 30 do
      world:step(1 / 60)
    end

    expect(world:getBody(weight).position.Y < 3).toBe(true)
  end)

  suite.it('lets jointed bodies overlap unless the joint says they collide', function()
    local function contacts(collide: boolean): number
      local world = physik.newWorld({ gravity = vector.create(0, 0, 0) })
      local left = world:addBody({ shape = physik.sphere(0.5) })
      local right = world:addBody({ shape = physik.sphere(0.5), position = vector.create(0.6, 0, 0) })
      world:addBallJoint(left, right, { collideConnected = collide })

      local hinge = world:addHingeJoint(left, right, { collideConnected = collide, motorSpeed = 1, maxMotorTorque = 2 })
      local count = world:step(1 / 60)
      world:removeJoint(hinge)

      return count
    end

    expect(contacts(false)).toBe(0)
    expect(contacts(true)).toBe(1)
  end)

  suite.it('lets the two bodies of a world\'s only joint overlap until it is removed', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0) })
    local left = world:addBody({ shape = physik.sphere(0.5) })
    local right = world:addBody({ shape = physik.sphere(0.5), position = vector.create(0.6, 0, 0) })
    local joint = world:addDistanceJoint(left, right, { length = 0.6 })

    expect(world:step(1 / 60)).toBe(0)
    world:removeJoint(joint)
    expect(world:step(1 / 60)).toBe(1)
  end)

  suite.it('finds a sphere in the cell it has crossed into, along each axis', function()
    local target = vector.create(0.5, 0.5, 0.5)

    for _, heading in { vector.create(1, 0, 0), vector.create(0, 1, 0), vector.create(0, 0, 1) } do
      local world = physik.newWorld({ gravity = vector.create(0, 0, 0) })
      world:addBody({ type = 'static', shape = physik.sphere(0.3), position = target })

      local mover = world:addBody({
        shape = physik.sphere(0.3),
        position = target - heading * 4,
        velocity = heading * 3,
        linearDamping = 0,
      })

      local contacts = 0

      for _ = 1, 90 do
        contacts += world:step(1 / 60)
      end

      -- Four cells on, it has met the sphere in its way, and has not gone through it.
      expect(contacts > 0).toBe(true)
      expect(vector.dot(world:getBody(mover).position - target, heading) < 0).toBe(true)
    end
  end)

  suite.it('pairs spheres too wide for a cell with what is two cells away', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0) })
    world:addBody({ type = 'static', shape = physik.sphere(0.9), position = vector.create(0.5, 0.5, 0.5) })
    world:addBody({ shape = physik.sphere(0.9), position = vector.create(2.2, 0.5, 0.5) })

    expect(world:step(1 / 60)).toBe(1)
  end)

  suite.it('turns a piece of a compound by its own rotation, then by its body\'s', function()
    local quarter = math.pi / 2
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0) })

    -- A capsule lies along Y. The piece's quarter turn about X lays it along Z, and the body's
    -- about Z leaves it there; turned in the other order it would lie along X.
    world:addBody({
      type = 'static',
      rotation = physik.axisAngle(vector.create(0, 0, 1), quarter),
      shape = physik.compound({
        { shape = physik.capsule(0.2, 0.6), rotation = physik.axisAngle(vector.create(1, 0, 0), quarter) },
      }),
    })

    local atItsEnd = world:addBody({ shape = physik.sphere(0.1), position = vector.create(0, 0, 0.85) })
    world:addBody({ shape = physik.sphere(0.1), position = vector.create(0.85, 0, 0) })

    expect(world:step(1 / 60)).toBe(1)

    local contact = world:contacts()[1]

    expect(contact.b).toBe(atItsEnd)
    expect(contact.normal.Z > 0.9).toBe(true)
  end)

  suite.it('wakes a sleeping body that a kinematic one runs into', function()
    local world = physik.newWorld()
    ground(world)

    -- Added before what will push it, so the sleeper is the first body of their pair.
    local ball = world:addBody({ shape = physik.sphere(0.4), position = vector.create(0, 0.4, 0) })

    for _ = 1, 120 do
      world:step(1 / 60)
    end

    expect(world:getBody(ball).sleeping).toBe(true)

    world:addBody({
      type = 'kinematic',
      shape = physik.box(vector.create(0.2, 0.4, 0.4)),
      position = vector.create(-2, 0.45, 0),
      velocity = vector.create(2, 0, 0),
    })

    for _ = 1, 60 do
      world:step(1 / 60)
    end

    expect(world:getBody(ball).sleeping).toBe(false)
    expect(world:getBody(ball).position.X > 0.2).toBe(true)
  end)
end)

suite.describe('queries', function()
  --- One of everything a ray can hit, a row along X, each named by its user data.
  local function gallery()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0), broadphaseCellSize = 4 })

    world:addBody({ type = 'static', shape = physik.sphere(1), position = vector.create(0, 0, 0), userData = 'sphere' })
    world:addBody({ type = 'static', shape = physik.box(vector.create(1, 0.5, 1)), position = vector.create(10, 0, 0), userData = 'box' })
    world:addBody({ type = 'static', shape = physik.capsule(0.5, 1), position = vector.create(20, 0, 0), userData = 'capsule' })

    world:addBody({
      type = 'static',
      position = vector.create(30, 0, 0),
      userData = 'compound',
      shape = physik.compound({
        { shape = physik.box(vector.create(1, 0.25, 1)) },
        { shape = physik.sphere(0.5), position = vector.create(0, 1, 0) },
      }),
    })

    world:addBody({
      type = 'static',
      position = vector.create(40, 0, 0),
      userData = 'mesh',
      shape = physik.triangleMesh(
        { vector.create(-2, 0.5, -2), vector.create(2, 0.5, -2), vector.create(2, 0.5, 2), vector.create(-2, 0.5, 2) },
        { 1, 3, 2, 1, 4, 3 }
      ),
    })

    world:addBody({
      type = 'static',
      position = vector.create(50, 0, 0),
      userData = 'heightfield',
      shape = physik.heightfield(3, 3, table.create(9, 0.75), { cellSize = 2, origin = vector.create(-2, 0, -2) }),
    })

    return world
  end

  local tops = { sphere = 1, box = 0.5, capsule = 1.5, compound = 1.5, mesh = 0.5, heightfield = 0.75 }
  local columns = { sphere = 0, box = 10, capsule = 20, compound = 30, mesh = 40, heightfield = 50 }

  suite.it('finds where a ray first meets each kind of shape', function()
    local world = gallery()

    for name, column in columns do
      local hit = world:raycast(vector.create(column + 0.1, 8, 0.1), vector.create(0, -2, 0), 20)

      expect(hit.userData).toBe(name)
      near(hit.distance, 8 - tops[name], 0.03)
      near(hit.point.Y, tops[name], 0.03)
      expect(hit.normal.Y > 0.9).toBe(true)
    end

    expect(world:raycast(vector.create(5, 8, 0), vector.create(0, -1, 0), 20)).toBe(nil)
    expect(world:raycast(vector.create(0, 8, 0), vector.create(0, -1, 0), 6.5)).toBe(nil)
    expect(world:raycast(vector.create(0, 8, 0), vector.create(0, 0, 0), 20)).toBe(nil)
    expect(world:raycast(vector.create(0, 8, 0), vector.create(0, -1, 0), 20, 2)).toBe(nil)
  end)

  suite.it('finds the nearest of the bodies along a ray, bounded or not', function()
    local world = gallery()
    local along = world:raycast(vector.create(-5, 0, 0), vector.create(1, 0, 0))
    local back = world:raycast(vector.create(60, 0.2, 0), vector.create(-1, 0, 0), 100)

    expect(along.userData).toBe('sphere')
    near(along.distance, 4, 1e-4)
    nearVector(along.normal, vector.create(-1, 0, 0), 1e-4)
    -- The heightfield and the mesh are surfaces: a ray below them passes under.
    expect(back.userData).toBe('compound')
    near(back.distance, 29, 1e-4)
  end)

  suite.it('finds a sphere where it has moved to, not where it was added', function()
    for _, broadphase in { 'grid', 'sap' } do
      local world = physik.newWorld({ gravity = vector.create(0, 0, 0), broadphase = broadphase, broadphaseCellSize = 2 })

      world:addBody({
        shape = physik.sphere(0.5),
        velocity = vector.create(1, 0, 0),
        linearDamping = 0,
        userData = 'ball',
      })

      for _ = 1, 300 do
        world:step(1 / 60)
      end

      local hit = world:raycast(vector.create(5, 5, 0), vector.create(0, -1, 0), 10)
      local swept = world:sphereCast(vector.create(5, 5, 0), vector.create(5, -5, 0), 0.2)

      expect(hit.userData).toBe('ball')
      near(hit.distance, 4.5, 1e-3)
      expect(swept.userData).toBe('ball')
      expect(world:raycast(vector.create(0, 5, 0), vector.create(0, -1, 0), 10)).toBe(nil)
      expect(named(world, world:queryAabb(vector.create(4.8, -1, -1), vector.create(5.2, 1, 1)), 1)).toEqual({ 'ball' })
    end
  end)

  suite.it('finds nothing along a ray that starts where no number is', function()
    local nowhere = 0 / 0
    local origins = { vector.create(nowhere, 8, 0), vector.create(0, nowhere, 0), vector.create(0, 8, nowhere) }

    for _, broadphase in { 'grid', 'sap' } do
      local world = physik.newWorld({ gravity = vector.create(0, 0, 0), broadphase = broadphase, broadphaseCellSize = 4 })

      -- Spheres small enough to be in the grid's cells, which is where such a ray is not.
      world:addBody({ type = 'static', shape = physik.sphere(1) })
      world:addBody({ type = 'static', shape = physik.sphere(0.5), position = vector.create(3, -2, 3) })

      for _, origin in origins do
        expect(world:raycast(origin, vector.create(0, -1, 0), 20)).toBe(nil)
        expect(world:raycast(origin, vector.create(1, -1, 1), 20)).toBe(nil)
        expect(world:raycast(origin, vector.create(0, -1, 0))).toBe(nil)
      end

      local drained, count = world:raycastBatchSoA(origins, [for _ in origins => vector.create(0, -1, 0)], { 20, 20, 20 })

      expect(count).toBe(3)
      expect(drained.hit).toEqual({ false, false, false })
      expect(world:raycast(vector.create(0, 8, 0), vector.create(0, -1, 0), 20).distance).toBe(7)
    end
  end)

  suite.it('answers a batch of rays as records and as columns, like the rays one by one', function()
    local world = gallery()
    local origins = [for column in range(0, 52.5, 2.5) => vector.create(column, 6, 0.2)]
    local directions = [for _ in origins => vector.create(0.05, -1, 0)]
    local reaches = [for index in range(1, #origins) => 5 + index % 3]
    local records = world:raycastBatch(origins, directions, reaches)
    local drained, count = world:raycastBatchSoA(origins, directions, reaches)

    expect(count).toBe(#origins)
    expect(records.count).toBe(#origins)
    expect(drained.count).toBe(#origins)

    local hits = 0

    for index, origin in origins do
      local hit = world:raycast(origin, directions[index], reaches[index])

      expect(drained.hit[index]).toBe(hit ~= nil)
      expect(records[index]).toEqual(hit)
      if not hit then continue end

      hits += 1
      expect(drained.handles[index]).toBe(hit.handle)
      expect(drained.distances[index]).toBe(hit.distance)
      expect(drained.points[index]).toBe(hit.point)
      expect(drained.normals[index]).toBe(hit.normal)
      expect(drained.userData[index]).toBe(hit.userData)
    end

    expect(hits).toBe(5)

    -- A shorter batch into the same tables leaves nothing of the longer one behind.
    world:raycastBatch({ origins[2] }, { directions[2] }, nil, nil, records)
    world:raycastBatchSoA({ origins[2] }, { directions[2] }, nil, nil, drained)
    expect(records.count).toBe(1)
    expect(#records).toBe(0)
    expect(drained.count).toBe(1)
    expect(#drained.hit).toBe(1)
    expect(#drained.handles).toBe(0)
  end)

  suite.it('finds the bodies a box overlaps, through the grid and through a scan', function()
    local world = gallery()
    local few, fewCount = world:queryAabb(vector.create(8, -1, -1), vector.create(21, 1, 1))

    expect(named(world, few, fewCount)).toEqual({ 'box', 'capsule' })
    expect(#few).toBe(2)

    local every, everyCount = world:queryAabb(vector.create(-500, -500, -500), vector.create(500, 500, 500), nil, few)

    expect(every).toBe(few)
    expect(named(world, every, everyCount)).toEqual({ 'box', 'capsule', 'compound', 'heightfield', 'mesh', 'sphere' })

    local none, noneCount = world:queryAabb(vector.create(4, 4, 4), vector.create(5, 5, 5), nil, few)

    expect(noneCount).toBe(0)
    expect(#none).toBe(0)
    expect(select(2, world:queryAabb(vector.create(-1, -1, -1), vector.create(1, 1, 1), 2))).toBe(0)

    expect(function()
      world:queryAabb(0, 1)
    end).toThrow('physik.queryAabb: bounds must be vectors')
  end)

  suite.it('finds through the grid what testing every body finds, in a crowd and where cells stop packing', function()
    -- Half a unit to a cell ends the cells that pack 16,384 units out along an axis, so the
    -- second crowd straddles that: some of its bodies are in cells and the rest are not.
    for _, centre in { vector.create(0, 0, 0), vector.create(16382, 0, -16382) } do
      local options = { gravity = vector.create(0, 0, 0), broadphaseCellSize = 0.5 }
      local walked = physik.newWorld(options)

      -- A limit of one cell sends every ray and every box to the test of each body in turn.
      local scanned = physik.newWorld({
        gravity = options.gravity,
        broadphaseCellSize = options.broadphaseCellSize,
        queryRayCellLimit = 1,
        queryAabbCellLimit = 1,
      })

      local state = 7

      local function random(low: number, high: number): number
        state = (state * 1103515245 + 12345) % 2147483648

        return low + (high - low) * state / 2147483648
      end

      local function scattered(reach: number): vector
        return vector.create(random(-reach, reach), random(-reach, reach), random(-reach, reach))
      end

      for index = 1, 400 do
        local kind = index % 4
        local shape

        if kind == 0 then
          shape = physik.sphere(random(0.05, 0.24))
        elseif kind == 1 then
          shape = physik.box(vector.create(random(0.05, 0.16), random(0.05, 0.16), random(0.05, 0.16)))
        elseif kind == 2 then
          shape = physik.capsule(random(0.05, 0.12), random(0, 0.1))
        else
          -- One in four is too wide for a cell, and is kept beside the grid.
          shape = physik.box(vector.create(random(0.3, 0.6), random(0.1, 0.3), random(0.1, 0.3)))
        end

        local body = { type = 'static', shape = shape, position = centre + scattered(4), userData = index }

        walked:addBody(body)
        scanned:addBody(body)
      end

      local hits, overlaps = 0, 0

      for _ = 1, 300 do
        local origin, direction, reach = centre + scattered(5), scattered(1), random(1, 12)
        local found, wanted = walked:raycast(origin, direction, reach), scanned:raycast(origin, direction, reach)

        expect(found ~= nil).toBe(wanted ~= nil)

        if found and wanted then
          hits += 1
          expect(found.handle).toBe(wanted.handle)
          expect(found.distance).toBe(wanted.distance)
        end

        local half = vector.create(random(0.2, 1.5), random(0.2, 1.5), random(0.2, 1.5))
        local inCells, cellCount = walked:queryAabb(origin - half, origin + half)
        local inTurn, turnCount = scanned:queryAabb(origin - half, origin + half)

        table.sort(inCells)
        table.sort(inTurn)
        overlaps += cellCount
        expect(cellCount).toBe(turnCount)
        expect(inCells).toEqual(inTurn)
      end

      -- The rays and the boxes did find things, so agreeing on them says something.
      expect(hits > 30).toBe(true)
      expect(overlaps > 300).toBe(true)
    end
  end)

  suite.it('counts a sphere\'s margin when a box asks what it overlaps', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0), contactMargin = 0.1 })

    -- One sphere too wide for a cell and one in a cell, each reaching 0.1 past its radius.
    world:addBody({ type = 'static', shape = physik.sphere(1) })
    world:addBody({ type = 'static', shape = physik.sphere(0.25), position = vector.create(0, 5, 0) })

    expect(select(2, world:queryAabb(vector.create(1.05, -0.1, -0.1), vector.create(2, 0.1, 0.1)))).toBe(1)
    expect(select(2, world:queryAabb(vector.create(1.15, -0.1, -0.1), vector.create(2, 0.1, 0.1)))).toBe(0)
    expect(select(2, world:queryAabb(vector.create(0.3, 4.9, -0.1), vector.create(1, 5.1, 0.1)))).toBe(1)
    expect(select(2, world:queryAabb(vector.create(0.4, 4.9, -0.1), vector.create(1, 5.1, 0.1)))).toBe(0)
  end)

  suite.it('sweeps a sphere to a body in a slot that another body had before it', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0) })

    world:addBody({ type = 'static', shape = physik.sphere(0.5), position = vector.create(0, 10, 0), userData = 'aside' })
    world:removeBody(world:addBody({ type = 'static', shape = physik.sphere(0.5), userData = 'gone' }))

    local again = world:addBody({ type = 'static', shape = physik.box(vector.create(0.5, 0.5, 0.5)), userData = 'again' })
    local hit = world:sphereCast(vector.create(0, 5, 0), vector.create(0, -5, 0), 0.25)

    expect(again).toBe(2 * Generation + 2)
    expect(hit.handle).toBe(again)
    expect(hit.userData).toBe('again')
    near(hit.center.Y, 0.75, 1e-3)
    expect(world:sphereCast(vector.create(0, 5, 0), vector.create(0, -5, 0), 0.25, nil, again)).toBe(nil)
  end)

  suite.it('sweeps a sphere to where it first touches each kind of shape', function()
    local world = gallery()

    for name, column in columns do
      local hit = world:sphereCast(vector.create(column + 0.1, 8, 0.1), vector.create(column + 0.1, -2, 0.1), 0.25)

      expect(hit.userData).toBe(name)
      near(hit.center.Y, tops[name] + 0.25, 0.03)
      near(hit.point.Y, tops[name], 0.03)
      near(hit.distance, 8 - tops[name] - 0.25, 0.03)
      near(hit.fraction, hit.distance / 10, 1e-6)
      expect(hit.normal.Y > 0.9).toBe(true)
    end

    expect(world:sphereCast(vector.create(5, 8, 0), vector.create(5, -2, 0), 0.25)).toBe(nil)

    expect(function()
      world:sphereCast(vector.create(5, 8, 0), vector.create(5, -2, 0), 0)
    end).toThrow('physik.sphereCast: expected start/end vectors and radius > 0')
  end)

  suite.it('answers a batch of sphere casts as columns, like the casts one by one', function()
    local world = gallery()
    local starts = [for column in range(0, 52.5, 2.5) => vector.create(column, 6, 0.2)]
    local ends = [for start in starts => start - vector.create(0, 7, 0)]
    local radii = [for index in range(1, #starts) => 0.2 + index % 2 * 0.2]
    local sphere = world:queryAabb(vector.create(-1, -1, -1), vector.create(1, 1, 1))[1]
    local ignored = [for _ in starts => sphere]
    local drained, count = world:sphereCastBatchSoA(starts, ends, radii, nil, ignored)
    local hits = 0

    expect(count).toBe(#starts)

    for index, start in starts do
      local hit = world:sphereCast(start, ends[index], radii[index], nil, sphere)

      expect(drained.hit[index]).toBe(hit ~= nil)
      if not hit then continue end

      hits += 1
      expect(drained.handles[index]).toBe(hit.handle)
      expect(drained.fractions[index]).toBe(hit.fraction)
      expect(drained.distances[index]).toBe(hit.distance)
      expect(drained.centers[index]).toBe(hit.center)
      expect(drained.points[index]).toBe(hit.point)
      expect(drained.normals[index]).toBe(hit.normal)
      expect(drained.userData[index]).toBe(hit.userData)
    end

    expect(hits).toBe(5)

    world:sphereCastBatchSoA({ starts[2] }, { ends[2] }, 0.5, nil, nil, drained)
    expect(drained.count).toBe(1)
    expect(#drained.hit).toBe(1)
    expect(#drained.handles).toBe(0)
  end)
end)

suite.describe('drains', function()
  local function scene()
    local world = physik.newWorld({ contactEvents = true })
    local floor = ground(world)

    local resting = world:addBody({ shape = physik.sphere(0.5), position = vector.create(0, 0.5, 0), userData = 'resting' })
    local falling = world:addBody({ shape = physik.sphere(0.5), position = vector.create(4, 30, 0), userData = 'falling' })
    local carried = world:addBody({ type = 'kinematic', shape = physik.sphere(0.5), position = vector.create(-4, 3, 0), userData = 'carried' })

    for _ = 1, 90 do
      world:step(1 / 60)
    end

    return world, floor, resting, falling, carried
  end

  suite.it('writes the transforms it is asked for and clears what a longer drain left', function()
    local world, floor, resting, falling, carried = scene()
    local state, count = world:writeTransforms()

    expect(count).toBe(4)
    expect(state.count).toBe(4)
    expect(state.handles).toEqual({ floor, resting, falling, carried })
    expect(state.positions[4]).toBe(vector.create(-4, 3, 0))
    expect(state.rotationW[1]).toBe(1)
    expect(state.rotations[1]).toBe(vector.create(0, 0, 0))
    expect(state.velocities[3].Y < -5).toBe(true)
    expect(state.angularVelocities[1]).toBe(vector.create(0, 0, 0))

    local dynamic, dynamicCount = world:writeTransforms(state, { dynamicOnly = true })

    expect(dynamic).toBe(state)
    expect(dynamicCount).toBe(2)
    expect(state.handles).toEqual({ resting, falling })
    expect(#state.positions).toBe(2)
    expect(#state.rotationW).toBe(2)

    local awake, awakeCount = world:writeTransforms(state, { dynamicOnly = true, awakeOnly = true })

    expect(awakeCount).toBe(1)
    expect(awake.handles).toEqual({ falling })
    expect(#state.angularVelocities).toBe(1)
  end)

  suite.it('writes linear state with the columns it is asked for', function()
    local world, floor, resting, falling, carried = scene()
    local state, count = world:writeLinearState()

    expect(count).toBe(4)
    expect(state.handles).toEqual({ floor, resting, falling, carried })
    expect(state.userData).toBe(nil)
    expect(#state.positions).toBe(4)
    expect(#state.velocities).toBe(4)

    local named, namedCount = world:writeLinearState(state, { awakeOnly = true, dynamicOnly = true, userData = true })

    expect(named).toBe(state)
    expect(namedCount).toBe(1)
    expect(state.handles).toEqual({ falling })
    expect(state.userData).toEqual({ 'falling' })
    expect(#state.positions).toBe(1)
    expect(#state.velocities).toBe(1)
  end)

  suite.it('writes contacts as columns and as records that agree', function()
    local world, floor, resting = scene()
    local state, count = world:writeContacts()
    local records = world:contacts()

    expect(count).toBe(1)
    expect(#records).toBe(1)
    expect(state.a).toEqual({ floor })
    expect(state.b).toEqual({ resting })
    expect(state.sensor).toEqual({ false })
    nearVector(state.normals[1], vector.create(0, 1, 0), 1e-4)
    near(state.points[1].Y, 0, 0.02)

    local record = records[1]

    expect(record.a).toBe(floor)
    expect(record.b).toBe(resting)
    expect(record.normal).toBe(state.normals[1])
    expect(record.point).toBe(state.points[1])
    expect(record.penetration).toBe(state.penetration[1])
    expect(record.normalImpulse).toBe(state.normalImpulse[1])
    expect(record.sensor).toBe(false)

    -- A drain that held more contacts is cut back to the ones there are.
    table.insert(state.a, 0)
    table.insert(state.sensor, true)
    world:writeContacts(state)
    expect(state.a).toEqual({ floor })
    expect(state.sensor).toEqual({ false })
    expect(state.count).toBe(1)
  end)

  suite.it('applies transforms back and moves the origin under every body', function()
    local world, floor, resting, falling, carried = scene()
    local state = world:writeTransforms()

    state.positions[2] = vector.create(7, 0.5, 7)
    state.velocities[3] = vector.create(0, 1, 0)
    expect(world:applyTransforms(state)).toBe(4)
    expect(world:getBody(resting).position).toBe(vector.create(7, 0.5, 7))
    expect(world:getBody(resting).sleeping).toBe(false)
    expect(world:getBody(falling).velocity).toBe(vector.create(0, 1, 0))
    expect(named(world, world:queryAabb(vector.create(6, 0, 6), vector.create(8, 1, 8)), 2)).toEqual({ 'ground', 'resting' })

    expect(world:shiftOrigin(vector.create(100, 0, -50))).toBe(world)
    expect(world:getBody(resting).position).toBe(vector.create(-93, 0.5, 57))
    expect(world:getBody(carried).position).toBe(vector.create(-104, 3, 50))
    expect(world:getBody(floor).position).toBe(vector.create(-100, -0.5, 50))
    expect(world:raycast(vector.create(-93, 5, 57), vector.create(0, -1, 0), 10).userData).toBe('resting')

    expect(function()
      world:applyTransforms({})
    end).toThrow('physik.applyTransforms: expected writeTransforms-compatible data')
  end)

  suite.it('moves the origin under the contacts of the last step too', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0) })
    world:addBody({ type = 'static', shape = physik.sphere(0.5) })
    world:addBody({ shape = physik.sphere(0.5), position = vector.create(0.9, 0, 0) })

    expect(world:step(1 / 60)).toBe(1)

    local shift = vector.create(8, -4, 2)
    local point = world:contacts()[1].point

    world:shiftOrigin(shift)
    expect(world:contacts()[1].point).toBe(point - shift)
    expect(world:writeContacts().points[1]).toBe(point - shift)
  end)

  suite.it('counts what the last step did', function()
    local world = scene()
    local stats = world:stats()

    expect(stats.frame).toBe(90)
    expect(stats.bodies).toBe(4)
    expect(stats.sleeping).toBe(1)
    expect(stats.contacts).toBe(1)
    expect(stats.sensorContacts).toBe(0)
    expect(stats.activeContacts).toBe(0)
    expect(stats.candidates).toBe(1)
    expect(stats.manifolds).toBe(1)
    expect(stats.joints).toBe(0)
    expect(stats.ccdHits).toBe(0)
    expect(stats.islands).toBe(2)
  end)
end)

return suite
