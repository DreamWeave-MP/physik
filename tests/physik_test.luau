local physik = require('@dream/physik')
local test = require('@jess/test')
local expect = test.expect
local suite = test.new()

suite.describe('linear drains', function()
  suite.it('emits only requested linear state and can identify rows by user data', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0) })
    world:addBody({
      shape = physik.sphere(1),
      position = vector.create(1, 2, 3),
      velocity = vector.create(4, 5, 6),
      mass = 1,
      userData = 77,
    })
    local state, count = world:writeLinearState({}, { dynamicOnly = true, userData = true, handles = false })
    expect(count).toBe(1)
    expect(state.handles).toBe(nil)
    expect(state.userData[1]).toBe(77)
    expect(state.positions[1].X).toBe(1)
    expect(state.velocities[1].Z).toBe(6)
  end)
end)

suite.describe('numeric grid broadphase', function()
  suite.it('finds contacts across adjacent occupied cells', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0), broadphase = 'grid', broadphaseCellSize = 1 })
    world:addBody({ shape = physik.sphere(0.4), position = vector.create(0.75, 0, 0), mass = 1 })
    world:addBody({ shape = physik.sphere(0.4), position = vector.create(1.25, 0, 0), mass = 1 })
    world:step(1 / 60)
    expect(world:stats().contacts > 0).toBe(true)
  end)

  suite.it('preserves the boundary fallback for cells at the packed-coordinate edge', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0), broadphase = 'grid', broadphaseCellSize = 1 })
    world:addBody({ shape = physik.sphere(0.4), position = vector.create(32766.75, 0, 0), mass = 1 })
    world:addBody({ shape = physik.sphere(0.4), position = vector.create(32767.25, 0, 0), mass = 1 })
    world:step(1 / 60)
    expect(world:stats().contacts > 0).toBe(true)
  end)
end)


suite.describe('coherent grid cache', function()
  suite.it('rebuilds occupied-neighbor topology when moving bodies cross cells', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0), broadphase = 'grid', broadphaseCellSize = 1 })
    world:addBody({ shape = physik.sphere(0.35), position = vector.create(0.2, 0, 0), velocity = vector.create(1.5, 0, 0), mass = 1, linearDamping = 0 })
    world:addBody({ shape = physik.sphere(0.35), position = vector.create(3.8, 0, 0), velocity = vector.create(-1.5, 0, 0), mass = 1, linearDamping = 0 })
    local touched = false
    for _ = 1, 120 do
      world:step(1 / 60)
      if world:stats().contacts > 0 then touched = true end
    end
    expect(touched).toBe(true)
  end)

  suite.it('keeps sphere AABB queries correct after coherent motion', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0), broadphase = 'grid', broadphaseCellSize = 2 })
    world:addBody({ shape = physik.sphere(0.5), position = vector.create(0, 0, 0), velocity = vector.create(1, 0, 0), mass = 1, linearDamping = 0 })
    for _ = 1, 30 do world:step(1 / 60) end
    local hits, count = world:queryAabb(vector.create(0.3, -1, -1), vector.create(1.2, 1, 1))
    expect(count).toBe(1)
    expect(#hits).toBe(1)
  end)
end)

suite.describe('linear state view', function()
  suite.it('exposes the current Physik-owned hot arrays without copying', function()
    local world = physik.newWorld({ gravity = vector.create(0, 0, 0) })
    world:addBody({ shape = physik.sphere(0.5), position = vector.create(1, 0, 0), velocity = vector.create(2, 0, 0), mass = 1, linearDamping = 0, userData = 9 })
    local view = world:linearStateView()
    expect(view.high).toBe(1)
    expect(view.userData[1]).toBe(9)
    expect(view.positions[1].X).toBe(1)
    world:step(0.5)
    expect(view.positions[1].X).toBe(2)
  end)
end)

return suite
