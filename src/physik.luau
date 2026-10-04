--!native
-- physik: compact 3D rigid-body physics for Jess/Luau.
-- Hot state uses Luau's native vector value type. Cold metadata stays scalar/packed-friendly.

local M = {}

local ZERO = vector.create(0, 0, 0)
local XAXIS = vector.create(1, 0, 0)
local YAXIS = vector.create(0, 1, 0)
local ZAXIS = vector.create(0, 0, 1)
local EPS = 1e-8
local HANDLE_STRIDE = 1048576 -- 20 index bits worth of exact-number handle space.

local BODY_STATIC = 0
local BODY_DYNAMIC = 1
local BODY_KINEMATIC = 2
local SHAPE_SPHERE = 1
local SHAPE_BOX = 2
local SHAPE_CAPSULE = 3

local CELL_BASE = 65536
local CELL_OFFSET = 32768
local CELL_BASE2 = 4294967296
local GRID_DX = {0, 1,1,1,1,1,1,1,1,1, 0,0,0, 0}
local GRID_DY = {0,-1,-1,-1,0,0,0,1,1,1, 1,1,1, 0}
local GRID_DZ = {0,-1, 0, 1,-1,0,1,-1,0,1,-1,0,1, 1}

local function packCell(x,y,z)
    if x < -CELL_OFFSET or x >= CELL_OFFSET or y < -CELL_OFFSET or y >= CELL_OFFSET or z < -CELL_OFFSET or z >= CELL_OFFSET then
        return nil
    end
    return (x + CELL_OFFSET) * CELL_BASE2 + (y + CELL_OFFSET) * CELL_BASE + (z + CELL_OFFSET)
end

M.STATIC = "static"
M.DYNAMIC = "dynamic"
M.KINEMATIC = "kinematic"

local function clamp(x, lo, hi)
    if x < lo then return lo end
    if x > hi then return hi end
    return x
end

local function signnz(x)
    if x < 0 then return -1 end
    return 1
end

local function safeUnit(v, fallback)
    local l2 = vector.dot(v, v)
    if l2 <= EPS then return fallback end
    return v / math.sqrt(l2)
end

local function qNormalize(qv, qw)
    local n2 = vector.dot(qv, qv) + qw * qw
    if n2 <= EPS then return ZERO, 1 end
    local inv = 1 / math.sqrt(n2)
    return qv * inv, qw * inv
end

local function qRotate(qv, qw, v)
    local t = vector.cross(qv, v) * 2
    return v + t * qw + vector.cross(qv, t)
end

local function qInvRotate(qv, qw, v)
    local nq = -qv
    local t = vector.cross(nq, v) * 2
    return v + t * qw + vector.cross(nq, t)
end

local function qAxes(qv, qw)
    local x, y, z = qv.X, qv.Y, qv.Z
    local xx, yy, zz = x*x, y*y, z*z
    local xy, xz, yz = x*y, x*z, y*z
    local wx, wy, wz = qw*x, qw*y, qw*z
    return
        vector.create(1 - 2*(yy + zz), 2*(xy + wz), 2*(xz - wy)),
        vector.create(2*(xy - wz), 1 - 2*(xx + zz), 2*(yz + wx)),
        vector.create(2*(xz + wy), 2*(yz - wx), 1 - 2*(xx + yy))
end

local function qIntegrate(qv, qw, omega, dt)
    local half = 0.5 * dt
    local nv = qv + (omega * qw + vector.cross(omega, qv)) * half
    local nw = qw - vector.dot(omega, qv) * half
    return qNormalize(nv, nw)
end

local function makeHandle(i, generation)
    return generation * HANDLE_STRIDE + i
end

local function splitHandle(handle)
    local i = handle % HANDLE_STRIDE
    return i, (handle - i) / HANDLE_STRIDE
end

local function bodyTypeCode(kind)
    if kind == nil or kind == "dynamic" then return BODY_DYNAMIC end
    if kind == "static" then return BODY_STATIC end
    if kind == "kinematic" then return BODY_KINEMATIC end
    error("physik: body type must be 'dynamic', 'static', or 'kinematic'", 3)
end

function M.sphere(radius)
    assert(type(radius) == "number" and radius > 0, "physik.sphere: radius must be > 0")
    return table.freeze({ kind = SHAPE_SPHERE, radius = radius })
end

function M.box(halfExtents)
    assert(typeof(halfExtents) == "vector", "physik.box: halfExtents must be a vector")
    assert(halfExtents.X > 0 and halfExtents.Y > 0 and halfExtents.Z > 0, "physik.box: extents must be > 0")
    return table.freeze({ kind = SHAPE_BOX, halfExtents = halfExtents })
end

function M.capsule(radius, halfHeight)
    assert(type(radius) == "number" and radius > 0, "physik.capsule: radius must be > 0")
    assert(type(halfHeight) == "number" and halfHeight >= 0, "physik.capsule: halfHeight must be >= 0")
    return table.freeze({ kind = SHAPE_CAPSULE, radius = radius, halfHeight = halfHeight })
end

function M.quat(x, y, z, w)
    local qv, qw = qNormalize(vector.create(x, y, z), w)
    return table.freeze({ v = qv, w = qw })
end

function M.axisAngle(axis, radians)
    axis = safeUnit(axis, YAXIS)
    local h = radians * 0.5
    return M.quat(axis.X * math.sin(h), axis.Y * math.sin(h), axis.Z * math.sin(h), math.cos(h))
end

local function inertiaFor(shape, mass)
    if mass <= 0 then return ZERO end
    if shape.kind == SHAPE_SPHERE then
        local i = 0.4 * mass * shape.radius * shape.radius
        local inv = if i > EPS then 1 / i else 0
        return vector.create(inv, inv, inv)
    elseif shape.kind == SHAPE_BOX then
        local h = shape.halfExtents
        local ix = (mass / 3) * (h.Y*h.Y + h.Z*h.Z)
        local iy = (mass / 3) * (h.X*h.X + h.Z*h.Z)
        local iz = (mass / 3) * (h.X*h.X + h.Y*h.Y)
        return vector.create(if ix > EPS then 1/ix else 0, if iy > EPS then 1/iy else 0, if iz > EPS then 1/iz else 0)
    else
        -- Capsule: compact approximation by a box enclosing the swept sphere.
        local r, hh = shape.radius, shape.halfHeight
        local hx, hy, hz = r, hh + r, r
        local ix = (mass / 3) * (hy*hy + hz*hz)
        local iy = (mass / 3) * (hx*hx + hz*hz)
        local iz = (mass / 3) * (hx*hx + hy*hy)
        return vector.create(if ix > EPS then 1/ix else 0, if iy > EPS then 1/iy else 0, if iz > EPS then 1/iz else 0)
    end
end

local World = {}
World.__index = World

local function applyInvInertia(self,i,v)
    if self.invMass[i]==0 then return ZERO end
    local ii=self.invInertia[i]
    -- Spheres are isotropic, so orientation cancels completely.
    if math.abs(ii.X-ii.Y)<1e-12 and math.abs(ii.X-ii.Z)<1e-12 then return v*ii.X end
    local a0,a1,a2=self.axis0[i],self.axis1[i],self.axis2[i]
    if a0==nil then
        a0,a1,a2=qAxes(self.qv[i],self.qw[i])
        self.axis0[i],self.axis1[i],self.axis2[i]=a0,a1,a2
    end
    return a0*(vector.dot(v,a0)*ii.X)+a1*(vector.dot(v,a1)*ii.Y)+a2*(vector.dot(v,a2)*ii.Z)
end

local function supportBox(center, qv, qw, half, dir)
    local a0, a1, a2 = qAxes(qv, qw)
    return center
        + a0 * (if vector.dot(dir, a0) >= 0 then half.X else -half.X)
        + a1 * (if vector.dot(dir, a1) >= 0 then half.Y else -half.Y)
        + a2 * (if vector.dot(dir, a2) >= 0 then half.Z else -half.Z)
end

local function supportBoxAxes(center,a0,a1,a2,half,dir)
    return center
        + a0 * (if vector.dot(dir,a0)>=0 then half.X else -half.X)
        + a1 * (if vector.dot(dir,a1)>=0 then half.Y else -half.Y)
        + a2 * (if vector.dot(dir,a2)>=0 then half.Z else -half.Z)
end

local function testBoxAxis(axis,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis)
    local l2=vector.dot(axis,axis)
    if l2<=1e-10 then return true,bestPen,bestAxis end
    local n=if math.abs(l2-1)<1e-5 then axis else axis/math.sqrt(l2)
    local ra=ah.X*math.abs(vector.dot(n,a0))+ah.Y*math.abs(vector.dot(n,a1))+ah.Z*math.abs(vector.dot(n,a2))
    local rb=bh.X*math.abs(vector.dot(n,b0))+bh.Y*math.abs(vector.dot(n,b1))+bh.Z*math.abs(vector.dot(n,b2))
    local dn=vector.dot(d,n)
    local pen=ra+rb-math.abs(dn)
    if pen<=0 then return false,bestPen,bestAxis end
    if pen<bestPen then
        bestPen=pen
        bestAxis=if dn>=0 then n else -n
    end
    return true,bestPen,bestAxis
end

local function closestPointSegment(p, a, b)
    local ab = b - a
    local d = vector.dot(ab, ab)
    if d <= EPS then return a end
    local t = clamp(vector.dot(p - a, ab) / d, 0, 1)
    return a + ab * t
end

local function closestSegmentSegment(p1, q1, p2, q2)
    local d1 = q1 - p1
    local d2 = q2 - p2
    local r = p1 - p2
    local a = vector.dot(d1, d1)
    local e = vector.dot(d2, d2)
    local f = vector.dot(d2, r)
    local s, t
    if a <= EPS and e <= EPS then
        return p1, p2
    elseif a <= EPS then
        s = 0
        t = clamp(f / e, 0, 1)
    else
        local c = vector.dot(d1, r)
        if e <= EPS then
            t = 0
            s = clamp(-c / a, 0, 1)
        else
            local b = vector.dot(d1, d2)
            local denom = a*e - b*b
            if math.abs(denom) > EPS then s = clamp((b*f - c*e) / denom, 0, 1) else s = 0 end
            t = (b*s + f) / e
            if t < 0 then
                t = 0
                s = clamp(-c / a, 0, 1)
            elseif t > 1 then
                t = 1
                s = clamp((b - c) / a, 0, 1)
            end
        end
    end
    return p1 + d1*s, p2 + d2*t
end

local function sphereSphere(pa, ra, pb, rb)
    local d = pb - pa
    local rr = ra + rb
    local d2 = vector.dot(d, d)
    if d2 >= rr*rr then return false end
    local dist = math.sqrt(math.max(d2, EPS))
    local n = if d2 > EPS then d / dist else YAXIS
    local penetration = rr - dist
    local point = pa + n * (ra - penetration * 0.5)
    return true, n, point, penetration
end

local function sphereBox(sp, radius, bp, bqv, bqw, half)
    local rel = qInvRotate(bqv, bqw, sp - bp)
    local closest = vector.create(
        clamp(rel.X, -half.X, half.X),
        clamp(rel.Y, -half.Y, half.Y),
        clamp(rel.Z, -half.Z, half.Z)
    )
    local d = rel - closest -- box -> sphere in box-local space.
    local d2 = vector.dot(d, d)
    if d2 > radius*radius then return false end
    local nLocal, penetration, pointLocal
    if d2 > EPS then
        local dist = math.sqrt(d2)
        nLocal = -(d / dist) -- sphere -> box
        penetration = radius - dist
        pointLocal = closest
    else
        local dx = half.X - math.abs(rel.X)
        local dy = half.Y - math.abs(rel.Y)
        local dz = half.Z - math.abs(rel.Z)
        if dx <= dy and dx <= dz then
            local outward = vector.create(signnz(rel.X), 0, 0)
            nLocal = -outward
            penetration = radius + dx
            pointLocal = rel + outward * dx
        elseif dy <= dz then
            local outward = vector.create(0, signnz(rel.Y), 0)
            nLocal = -outward
            penetration = radius + dy
            pointLocal = rel + outward * dy
        else
            local outward = vector.create(0, 0, signnz(rel.Z))
            nLocal = -outward
            penetration = radius + dz
            pointLocal = rel + outward * dz
        end
    end
    return true, qRotate(bqv, bqw, nLocal), bp + qRotate(bqv, bqw, pointLocal), penetration
end

local function boxBox(ap,a0,a1,a2,ah,bp,b0,b1,b2,bh)
    local d=bp-ap
    local bestPen,bestAxis=math.huge,nil
    local ok

    ok,bestPen,bestAxis=testBoxAxis(a0,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(a1,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(a2,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(b0,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(b1,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(b2,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(vector.cross(a0,b0),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(vector.cross(a0,b1),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(vector.cross(a0,b2),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(vector.cross(a1,b0),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(vector.cross(a1,b1),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(vector.cross(a1,b2),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(vector.cross(a2,b0),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(vector.cross(a2,b1),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end
    ok,bestPen,bestAxis=testBoxAxis(vector.cross(a2,b2),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis); if not ok then return false end

    local n=bestAxis or safeUnit(d,YAXIS)
    local sa=supportBoxAxes(ap,a0,a1,a2,ah,n)
    local sb=supportBoxAxes(bp,b0,b1,b2,bh,-n)
    return true,n,(sa+sb)*0.5,bestPen
end

local function capsuleEndpoints(p, qv, qw, hh)
    local axis = qRotate(qv, qw, YAXIS)
    local d = axis * hh
    return p - d, p + d
end

local function capsuleSphere(cp, cqv, cqw, cr, chh, sp, sr)
    local a,b = capsuleEndpoints(cp,cqv,cqw,chh)
    local q = closestPointSegment(sp,a,b)
    return sphereSphere(q,cr,sp,sr)
end

local function capsuleCapsule(ap,aqv,aqw,ar,ahh,bp,bqv,bqw,br,bhh)
    local a0,a1 = capsuleEndpoints(ap,aqv,aqw,ahh)
    local b0,b1 = capsuleEndpoints(bp,bqv,bqw,bhh)
    local ca,cb = closestSegmentSegment(a0,a1,b0,b1)
    return sphereSphere(ca,ar,cb,br)
end

local function capsuleBox(cp,cqv,cqw,cr,chh,bp,bqv,bqw,bh)
    -- Three sphere probes over the capsule spine. Conservative enough for character/controller-ish use,
    -- while keeping the hot narrow phase branch small. A segment/OBB closest-point kernel is a v0.2 target.
    local a,b = capsuleEndpoints(cp,cqv,cqw,chh)
    local m = (a+b)*0.5
    local ok,n,p,pen = sphereBox(a,cr,bp,bqv,bqw,bh)
    local bestOk,bestN,bestP,bestPen = ok,n,p,pen
    ok,n,p,pen = sphereBox(m,cr,bp,bqv,bqw,bh)
    if ok and (not bestOk or pen > bestPen) then bestOk,bestN,bestP,bestPen=ok,n,p,pen end
    ok,n,p,pen = sphereBox(b,cr,bp,bqv,bqw,bh)
    if ok and (not bestOk or pen > bestPen) then bestOk,bestN,bestP,bestPen=ok,n,p,pen end
    if not bestOk then return false end
    return true,bestN,bestP,bestPen
end

function World:_check(handle)
    assert(type(handle) == "number", "physik: body handle must be a number")
    local i,g = splitHandle(handle)
    assert(i >= 1 and self.alive[i] and self.generation[i] == g, "physik: stale or invalid body handle")
    return i
end

function World:_allocSlot()
    local n = #self.free
    local i
    if n > 0 then
        i = self.free[n]
        self.free[n] = nil
        self.generation[i] += 1
    else
        i = self.high + 1
        assert(i < HANDLE_STRIDE, "physik: world body capacity exceeded")
        self.high = i
        self.generation[i] = 1
    end
    self.alive[i] = true
    self.aliveCount += 1
    self.order[#self.order+1] = i
    return i
end

function World:addBody(desc)
    desc = desc or {}
    local shape = desc.shape or M.sphere(0.5)
    assert(type(shape) == "table" and (shape.kind == SHAPE_SPHERE or shape.kind == SHAPE_BOX or shape.kind == SHAPE_CAPSULE), "physik.addBody: invalid shape")
    local kind = bodyTypeCode(desc.type)
    local mass = if kind == BODY_DYNAMIC then (desc.mass or 1) else 0
    assert(mass >= 0, "physik.addBody: mass must be >= 0")
    local i = self:_allocSlot()
    self.bodyType[i] = kind
    self.shape[i] = shape
    self.pos[i] = desc.position or ZERO
    self.linVel[i] = desc.velocity or ZERO
    self.angVel[i] = desc.angularVelocity or ZERO
    self.force[i] = ZERO
    self.torque[i] = ZERO
    local r = desc.rotation
    if r == nil then
        self.qv[i], self.qw[i] = ZERO, 1
    elseif type(r) == "table" and typeof(r.v) == "vector" and type(r.w) == "number" then
        self.qv[i], self.qw[i] = qNormalize(r.v, r.w)
    else
        error("physik.addBody: rotation must come from physik.quat/axisAngle", 2)
    end
    self.axis0[i],self.axis1[i],self.axis2[i]=qAxes(self.qv[i],self.qw[i])
    self.invMass[i] = if kind == BODY_DYNAMIC and mass > EPS then 1/mass else 0
    self.invInertia[i] = if kind == BODY_DYNAMIC then inertiaFor(shape,mass) else ZERO
    self.friction[i] = desc.friction or 0.6
    self.restitution[i] = desc.restitution or 0.05
    self.linearDamping[i] = desc.linearDamping or 0.02
    self.angularDamping[i] = desc.angularDamping or 0.05
    self.group[i] = desc.group or 1
    self.mask[i] = if desc.mask == nil then 0xffffffff else desc.mask
    self.userData[i] = desc.userData
    self.sleeping[i] = false
    self.sleepTime[i] = 0
    self.aabbMin[i] = ZERO
    self.aabbMax[i] = ZERO
    return makeHandle(i,self.generation[i])
end

function World:removeBody(handle)
    local i = self:_check(handle)
    self.alive[i] = false
    self.aliveCount -= 1
    self.free[#self.free+1] = i
    self.userData[i] = nil
    self.sleeping[i] = false
end

function World:getBody(handle)
    local i = self:_check(handle)
    return {
        handle = handle,
        type = if self.bodyType[i] == BODY_STATIC then "static" elseif self.bodyType[i] == BODY_KINEMATIC then "kinematic" else "dynamic",
        position = self.pos[i],
        rotation = { v = self.qv[i], w = self.qw[i] },
        velocity = self.linVel[i],
        angularVelocity = self.angVel[i],
        sleeping = self.sleeping[i],
        userData = self.userData[i],
    }
end

function World:setPosition(handle, p, wake)
    local i = self:_check(handle)
    self.pos[i] = p
    if wake ~= false then self:wake(handle) end
end

function World:setRotation(handle,q,wake)
    local i=self:_check(handle)
    self.qv[i],self.qw[i]=qNormalize(q.v,q.w)
    self.axis0[i],self.axis1[i],self.axis2[i]=qAxes(self.qv[i],self.qw[i])
    if wake~=false then self:wake(handle) end
end

function World:setVelocity(handle, v)
    local i = self:_check(handle)
    self.linVel[i] = v
    self.sleeping[i] = false
    self.sleepTime[i] = 0
end

function World:setAngularVelocity(handle, v)
    local i = self:_check(handle)
    self.angVel[i] = v
    self.sleeping[i] = false
    self.sleepTime[i] = 0
end

function World:wake(handle)
    local i = self:_check(handle)
    self.sleeping[i] = false
    self.sleepTime[i] = 0
end

function World:applyForce(handle, f)
    local i = self:_check(handle)
    if self.bodyType[i] ~= BODY_DYNAMIC then return end
    self.force[i] += f
    self.sleeping[i] = false
    self.sleepTime[i] = 0
end

function World:applyForceAtPoint(handle, f, point)
    local i = self:_check(handle)
    if self.bodyType[i] ~= BODY_DYNAMIC then return end
    self.force[i] += f
    self.torque[i] += vector.cross(point - self.pos[i], f)
    self.sleeping[i] = false
    self.sleepTime[i] = 0
end

function World:applyImpulse(handle, impulse, point)
    local i = self:_check(handle)
    if self.bodyType[i] ~= BODY_DYNAMIC then return end
    self.linVel[i] += impulse * self.invMass[i]
    if point ~= nil then
        self.angVel[i] += applyInvInertia(self,i,vector.cross(point-self.pos[i],impulse))
    end
    self.sleeping[i] = false
    self.sleepTime[i] = 0
end

function World:addDistanceJoint(aHandle,bHandle,opts)
    local a,b = self:_check(aHandle),self:_check(bHandle)
    opts = opts or {}
    local la = opts.localAnchorA or ZERO
    local lb = opts.localAnchorB or ZERO
    local wa = self.pos[a] + qRotate(self.qv[a],self.qw[a],la)
    local wb = self.pos[b] + qRotate(self.qv[b],self.qw[b],lb)
    local id = self.nextJoint
    self.nextJoint += 1
    self.joints[id] = {
        a=a,b=b, ga=self.generation[a],gb=self.generation[b],
        la=la,lb=lb,
        rest=opts.length or vector.magnitude(wb-wa),
        stiffness=opts.stiffness or 0.25,
        damping=opts.damping or 0.05,
        impulse=0,
    }
    return id
end

function World:removeJoint(id)
    self.joints[id] = nil
end

function World:_updateAabb(i)
    local shape = self.shape[i]
    local p = self.pos[i]
    if shape.kind == SHAPE_SPHERE then
        local e = vector.create(shape.radius,shape.radius,shape.radius)
        self.aabbMin[i],self.aabbMax[i] = p-e,p+e
    elseif shape.kind == SHAPE_BOX then
        local a0,a1,a2=qAxes(self.qv[i],self.qw[i])
        self.axis0[i],self.axis1[i],self.axis2[i]=a0,a1,a2
        local h = shape.halfExtents
        local e = vector.create(
            math.abs(a0.X)*h.X + math.abs(a1.X)*h.Y + math.abs(a2.X)*h.Z,
            math.abs(a0.Y)*h.X + math.abs(a1.Y)*h.Y + math.abs(a2.Y)*h.Z,
            math.abs(a0.Z)*h.X + math.abs(a1.Z)*h.Y + math.abs(a2.Z)*h.Z
        )
        self.aabbMin[i],self.aabbMax[i] = p-e,p+e
    else
        local a0,axis,a2=qAxes(self.qv[i],self.qw[i])
        self.axis0[i],self.axis1[i],self.axis2[i]=a0,axis,a2
        local h = shape.halfHeight
        local r = shape.radius
        local e = vector.create(math.abs(axis.X)*h+r,math.abs(axis.Y)*h+r,math.abs(axis.Z)*h+r)
        self.aabbMin[i],self.aabbMax[i] = p-e,p+e
    end
end

function World:_sortOrder()
    -- Compact dead slots opportunistically, then insertion-sort on min X. Physics motion is coherent,
    -- so after the first frame this is normally close to O(n).
    local order = self.order
    if #order > self.aliveCount + 64 then
        local n = 0
        for k=1,#order do
            local i=order[k]
            if self.alive[i] then n+=1; order[n]=i end
        end
        for k=#order,n+1,-1 do order[k]=nil end
    end
    for k=2,#order do
        local v=order[k]
        local x=self.aabbMin[v].X
        local j=k-1
        while j>=1 and self.aabbMin[order[j]].X > x do
            order[j+1]=order[j]
            j-=1
        end
        order[j+1]=v
    end
end

local impulseDenom

function World:_addContact(a,b,n,p,pen)
    local c=self.contactCount+1
    self.contactCount=c
    self.ca[c],self.cb[c]=a,b
    self.cn[c],self.cp[c],self.cpen[c]=n,p,pen
    self.cnormalImpulse[c]=0
    self.ctangentImpulse1[c]=0
    self.ctangentImpulse2[c]=0
end

function World:_narrow(a,b)
    local sa,sb=self.shape[a],self.shape[b]
    local ap,bp=self.pos[a],self.pos[b]
    local aqv,aqw,bqv,bqw=self.qv[a],self.qw[a],self.qv[b],self.qw[b]
    local ok,n,p,pen
    if sa.kind==SHAPE_SPHERE and sb.kind==SHAPE_SPHERE then
        ok,n,p,pen=sphereSphere(ap,sa.radius,bp,sb.radius)
    elseif sa.kind==SHAPE_SPHERE and sb.kind==SHAPE_BOX then
        ok,n,p,pen=sphereBox(ap,sa.radius,bp,bqv,bqw,sb.halfExtents)
    elseif sa.kind==SHAPE_BOX and sb.kind==SHAPE_SPHERE then
        ok,n,p,pen=sphereBox(bp,sb.radius,ap,aqv,aqw,sa.halfExtents)
        if ok then n=-n end
    elseif sa.kind==SHAPE_BOX and sb.kind==SHAPE_BOX then
        ok,n,p,pen=boxBox(ap,self.axis0[a],self.axis1[a],self.axis2[a],sa.halfExtents,bp,self.axis0[b],self.axis1[b],self.axis2[b],sb.halfExtents)
    elseif sa.kind==SHAPE_CAPSULE and sb.kind==SHAPE_SPHERE then
        ok,n,p,pen=capsuleSphere(ap,aqv,aqw,sa.radius,sa.halfHeight,bp,sb.radius)
    elseif sa.kind==SHAPE_SPHERE and sb.kind==SHAPE_CAPSULE then
        ok,n,p,pen=capsuleSphere(bp,bqv,bqw,sb.radius,sb.halfHeight,ap,sa.radius)
        if ok then n=-n end
    elseif sa.kind==SHAPE_CAPSULE and sb.kind==SHAPE_CAPSULE then
        ok,n,p,pen=capsuleCapsule(ap,aqv,aqw,sa.radius,sa.halfHeight,bp,bqv,bqw,sb.radius,sb.halfHeight)
    elseif sa.kind==SHAPE_CAPSULE and sb.kind==SHAPE_BOX then
        ok,n,p,pen=capsuleBox(ap,aqv,aqw,sa.radius,sa.halfHeight,bp,bqv,bqw,sb.halfExtents)
    elseif sa.kind==SHAPE_BOX and sb.kind==SHAPE_CAPSULE then
        ok,n,p,pen=capsuleBox(bp,bqv,bqw,sb.radius,sb.halfHeight,ap,aqv,aqw,sa.halfExtents)
        if ok then n=-n end
    end
    if ok then self:_addContact(a,b,n,p,pen) end
end

function World:_tryPair(a,b)
    if a==b or not self.alive[a] or not self.alive[b] then return end
    if self.bodyType[a]==BODY_STATIC and self.bodyType[b]==BODY_STATIC then return end
    if bit32.band(self.mask[a],self.group[b])==0 or bit32.band(self.mask[b],self.group[a])==0 then return end
    local amin,amax=self.aabbMin[a],self.aabbMax[a]
    local bmin,bmax=self.aabbMin[b],self.aabbMax[b]
    if bmin.X>amax.X or bmax.X<amin.X or bmin.Y>amax.Y or bmax.Y<amin.Y or bmin.Z>amax.Z or bmax.Z<amin.Z then return end
    self.candidateCount+=1
    self:_narrow(a,b)
end

function World:_collideSap()
    self:_sortOrder()
    local order=self.order
    for oi=1,#order do
        local a=order[oi]
        if not self.alive[a] then continue end
        local amax=self.aabbMax[a]
        for oj=oi+1,#order do
            local b=order[oj]
            if not self.alive[b] then continue end
            if self.aabbMin[b].X > amax.X then break end
            self:_tryPair(a,b)
        end
    end
end

function World:_collideGrid()
    local heads,keys,nextIn=self.gridHeads,self.gridKeys,self.gridNext
    for k=1,#keys do heads[keys[k]]=nil; keys[k]=nil end
    local large=self.gridLarge
    for k=#large,1,-1 do large[k]=nil end
    local cs=self.broadphaseCellSize

    for i=1,self.high do
        if not self.alive[i] then continue end
        local mn,mx=self.aabbMin[i],self.aabbMax[i]
        local ex,ey,ez=mx.X-mn.X,mx.Y-mn.Y,mx.Z-mn.Z
        local p=self.pos[i]
        local cx,cy,cz=math.floor(p.X/cs),math.floor(p.Y/cs),math.floor(p.Z/cs)
        local key=if ex<=cs and ey<=cs and ez<=cs then packCell(cx,cy,cz) else nil
        if key then
            self.gridX[i],self.gridY[i],self.gridZ[i]=cx,cy,cz
            self.gridKey[i]=key
            local old=heads[key]
            if old==nil then
                local nk=#keys+1
                keys[nk]=key; self.gridKeyX[nk]=cx; self.gridKeyY[nk]=cy; self.gridKeyZ[nk]=cz
                old=0
            end
            nextIn[i]=old
            heads[key]=i
        else
            self.gridKey[i]=nil
            large[#large+1]=i
        end
    end

    -- Small-small: walk occupied cells, not bodies. Every cell checks only the lexicographic
    -- forward half of its 3x3x3 neighborhood, so pair generation is unique without a pair set.
    for ki=1,#keys do
        local key0=keys[ki]
        local cx,cy,cz=self.gridKeyX[ki],self.gridKeyY[ki],self.gridKeyZ[ki]
        local a=heads[key0]
        while a and a~=0 do
            local b=nextIn[a]
            while b and b~=0 do self:_tryPair(a,b); b=nextIn[b] end
            a=nextIn[a]
        end
        for o=2,14 do
            local key=packCell(cx+GRID_DX[o],cy+GRID_DY[o],cz+GRID_DZ[o])
            if key then
                a=heads[key0]
                while a and a~=0 do
                    local b2=heads[key]
                    while b2 and b2~=0 do self:_tryPair(a,b2); b2=nextIn[b2] end
                    a=nextIn[a]
                end
            end
        end
    end

    -- Bodies wider than one grid cell stay out of the hash grid. Pair anything involving
    -- them with a coherent X sweep instead of the old O(L^2) fallback. Small-small pairs
    -- are skipped here because the grid already emitted them.
    if #large>0 then
        self:_sortOrder()
        local order=self.order
        for oi=1,#order do
            local a=order[oi]
            if not self.alive[a] then continue end
            local amaxX=self.aabbMax[a].X
            local aLarge=self.gridKey[a]==nil
            for oj=oi+1,#order do
                local b=order[oj]
                if not self.alive[b] then continue end
                if self.aabbMin[b].X>amaxX then break end
                if aLarge or self.gridKey[b]==nil then self:_tryPair(a,b) end
            end
        end
    end
end

function World:_collide()
    self.contactCount=0
    self.candidateCount=0
    for i=1,self.high do if self.alive[i] then self:_updateAabb(i) end end
    if self.broadphase=="sap" then self:_collideSap() else self:_collideGrid() end
end

local function pointVelocity(self,i,r)
    if self.sleeping[i] then return ZERO end
    return self.linVel[i]+vector.cross(self.angVel[i],r)
end

impulseDenom = function(self,i,r,axis)
    if self.invMass[i]==0 or self.sleeping[i] then return 0 end
    local rx=vector.cross(r,axis)
    local ang=vector.cross(applyInvInertia(self,i,rx),r)
    return self.invMass[i]+vector.dot(axis,ang)
end

local function contactTangents(n)
    -- Pick the least-aligned cardinal axis to avoid degeneracy, then keep the
    -- basis fixed for the whole solve. This lets tangent effective mass be hoisted
    -- out of the velocity-iteration loop.
    local ax,ay,az=math.abs(n.X),math.abs(n.Y),math.abs(n.Z)
    local ref=if ax<=ay and ax<=az then XAXIS elseif ay<=az then YAXIS else ZAXIS
    local t1=safeUnit(vector.cross(n,ref),XAXIS)
    return t1,vector.cross(n,t1)
end

function World:_wakeContacts()
    local speed2=self.wakeSpeed*self.wakeSpeed
    local wakePen=self.wakePenetration
    for c=1,self.contactCount do
        local a,b=self.ca[c],self.cb[c]
        local sa,sb=self.sleeping[a],self.sleeping[b]
        if sa==sb then continue end
        local sleeper=if sa then a else b
        local mover=if sa then b else a
        if self.bodyType[sleeper]~=BODY_DYNAMIC then continue end
        local v=self.linVel[mover]
        local av=self.angVel[mover]
        if vector.dot(v,v)>speed2 or vector.dot(av,av)>speed2 or self.cpen[c]>wakePen then
            self.sleeping[sleeper]=false
            self.sleepTime[sleeper]=0
        end
    end
end

function World:_prepareContacts()
    for c=1,self.contactCount do
        local a,b=self.ca[c],self.cb[c]
        local p,n=self.cp[c],self.cn[c]
        local ra,rb=p-self.pos[a],p-self.pos[b]
        self.cra[c],self.crb[c]=ra,rb
        local d=impulseDenom(self,a,ra,n)+impulseDenom(self,b,rb,n)
        self.cnormalMass[c]=if d>EPS then 1/d else 0
        local t1,t2=contactTangents(n)
        self.ctangent1[c],self.ctangent2[c]=t1,t2
        local d1=impulseDenom(self,a,ra,t1)+impulseDenom(self,b,rb,t1)
        local d2=impulseDenom(self,a,ra,t2)+impulseDenom(self,b,rb,t2)
        self.ctangentMass1[c]=if d1>EPS then 1/d1 else 0
        self.ctangentMass2[c]=if d2>EPS then 1/d2 else 0
        self.cfriction[c]=math.sqrt(self.friction[a]*self.friction[b])
        self.crest[c]=math.max(self.restitution[a],self.restitution[b])
    end
end

local function applyPairImpulse(self,a,b,ra,rb,impulse)
    if self.invMass[a]>0 and not self.sleeping[a] then
        self.linVel[a]-=impulse*self.invMass[a]
        local spin=vector.cross(ra,impulse)
        if vector.dot(spin,spin)>EPS then self.angVel[a]-=applyInvInertia(self,a,spin) end
    end
    if self.invMass[b]>0 and not self.sleeping[b] then
        self.linVel[b]+=impulse*self.invMass[b]
        local spin=vector.cross(rb,impulse)
        if vector.dot(spin,spin)>EPS then self.angVel[b]+=applyInvInertia(self,b,spin) end
    end
end

function World:_solveContacts(dt)
    local iterations=self.velocityIterations
    for iter=1,iterations do
        for c=1,self.contactCount do
            local a,b=self.ca[c],self.cb[c]
            if not self.alive[a] or not self.alive[b] then continue end
            local n=self.cn[c]
            local ra,rb=self.cra[c],self.crb[c]
            local rv=pointVelocity(self,b,rb)-pointVelocity(self,a,ra)
            local vn=vector.dot(rv,n)
            local normalMass=self.cnormalMass[c]
            if normalMass>0 then
                local restitution=self.crest[c]
                local bounce=if iter==1 and vn < -self.restitutionThreshold then -restitution*vn else 0
                local lambda=(bounce-vn)*normalMass
                local old=self.cnormalImpulse[c]
                local now=math.max(old+lambda,0)
                lambda=now-old
                self.cnormalImpulse[c]=now
                if lambda~=0 then
                    applyPairImpulse(self,a,b,ra,rb,n*lambda)
                end
            end

            rv=pointVelocity(self,b,rb)-pointVelocity(self,a,ra)
            local t1,t2=self.ctangent1[c],self.ctangent2[c]
            local m1,m2=self.ctangentMass1[c],self.ctangentMass2[c]
            if m1>0 or m2>0 then
                local old1,old2=self.ctangentImpulse1[c],self.ctangentImpulse2[c]
                local new1=old1-vector.dot(rv,t1)*m1
                local new2=old2-vector.dot(rv,t2)*m2
                local limit=self.cfriction[c]*self.cnormalImpulse[c]
                local mag2=new1*new1+new2*new2
                if mag2>limit*limit and mag2>EPS then
                    local scale=limit/math.sqrt(mag2)
                    new1,new2=new1*scale,new2*scale
                end
                local d1,d2=new1-old1,new2-old2
                self.ctangentImpulse1[c],self.ctangentImpulse2[c]=new1,new2
                if d1~=0 or d2~=0 then applyPairImpulse(self,a,b,ra,rb,t1*d1+t2*d2) end
            end
        end
        self:_solveJoints(dt)
    end
end

function World:_solveJoints(dt)
    for id,j in self.joints do
        local a,b=j.a,j.b
        if not self.alive[a] or not self.alive[b] or self.generation[a]~=j.ga or self.generation[b]~=j.gb then
            self.joints[id]=nil
            continue
        end
        local ra=qRotate(self.qv[a],self.qw[a],j.la)
        local rb=qRotate(self.qv[b],self.qw[b],j.lb)
        local wa=self.pos[a]+ra
        local wb=self.pos[b]+rb
        local d=wb-wa
        local len=vector.magnitude(d)
        local n=if len>EPS then d/len else YAXIS
        local rv=pointVelocity(self,b,rb)-pointVelocity(self,a,ra)
        local denom=impulseDenom(self,a,ra,n)+impulseDenom(self,b,rb,n)
        if denom>EPS then
            local error=len-j.rest
            local bias=(error*j.stiffness)/math.max(dt,1e-5)
            local lambda=-(vector.dot(rv,n)+bias+j.damping*vector.dot(rv,n))/denom
            applyPairImpulse(self,a,b,ra,rb,n*lambda)
        end
    end
end

function World:_correctPositions()
    local slop=self.penetrationSlop
    local percent=self.positionCorrection
    for c=1,self.contactCount do
        local a,b=self.ca[c],self.cb[c]
        local imA=if self.sleeping[a] then 0 else self.invMass[a]
        local imB=if self.sleeping[b] then 0 else self.invMass[b]
        local sum=imA+imB
        if sum<=EPS then continue end
        local depth=math.max(self.cpen[c]-slop,0)
        if depth<=0 then continue end
        local correction=self.cn[c]*(depth*percent/sum)
        if imA>0 then self.pos[a]-=correction*imA end
        if imB>0 then self.pos[b]+=correction*imB end
    end
end

function World:_integrate(dt)
    local g=self.gravity
    for i=1,self.high do
        if not self.alive[i] then continue end
        local typ=self.bodyType[i]
        if typ==BODY_DYNAMIC then
            if self.sleeping[i] then
                self.force[i],self.torque[i]=ZERO,ZERO
                continue
            end
            local im=self.invMass[i]
            self.linVel[i]+=(g+self.force[i]*im)*dt
            local tq=self.torque[i]
            if vector.dot(tq,tq)>EPS then self.angVel[i]+=applyInvInertia(self,i,tq)*dt end
            self.linVel[i]*=1/(1+self.linearDamping[i]*dt)
            self.angVel[i]*=1/(1+self.angularDamping[i]*dt)
            self.pos[i]+=self.linVel[i]*dt
            local av=self.angVel[i]
            if vector.dot(av,av)>EPS then self.qv[i],self.qw[i]=qIntegrate(self.qv[i],self.qw[i],av,dt) end
            self.force[i],self.torque[i]=ZERO,ZERO
        elseif typ==BODY_KINEMATIC then
            self.pos[i]+=self.linVel[i]*dt
            local av=self.angVel[i]
            if vector.dot(av,av)>EPS then self.qv[i],self.qw[i]=qIntegrate(self.qv[i],self.qw[i],av,dt) end
        end
    end
end

function World:_sleep(dt)
    if not self.allowSleep then return end
    local lin2=self.sleepLinear*self.sleepLinear
    local ang2=self.sleepAngular*self.sleepAngular
    for i=1,self.high do
        if self.alive[i] and self.bodyType[i]==BODY_DYNAMIC and not self.sleeping[i] then
            if vector.dot(self.linVel[i],self.linVel[i])<lin2 and vector.dot(self.angVel[i],self.angVel[i])<ang2 then
                self.sleepTime[i]+=dt
                if self.sleepTime[i]>=self.sleepDelay then
                    self.sleeping[i]=true
                    self.linVel[i],self.angVel[i]=ZERO,ZERO
                end
            else
                self.sleepTime[i]=0
            end
        end
    end
end

function World:step(dt,substeps)
    assert(type(dt)=="number" and dt>0,"physik.step: dt must be > 0")
    substeps=substeps or 1
    assert(substeps>=1 and substeps%1==0,"physik.step: substeps must be a positive integer")
    local h=dt/substeps
    for _=1,substeps do
        self:_integrate(h)
        self:_collide()
        self:_wakeContacts()
        self:_prepareContacts()
        self:_solveContacts(h)
        self:_correctPositions()
        self:_sleep(h)
    end
    self.frame+=1
    return self.contactCount
end

local function raySphere(origin,dir,maxDist,center,radius)
    local m=origin-center
    local b=vector.dot(m,dir)
    local c=vector.dot(m,m)-radius*radius
    if c>0 and b>0 then return nil end
    local discr=b*b-c
    if discr<0 then return nil end
    local t=-b-math.sqrt(discr)
    if t<0 then t=0 end
    if t>maxDist then return nil end
    local p=origin+dir*t
    return t,safeUnit(p-center,-dir),p
end

local function rayBox(origin,dir,maxDist,center,qv,qw,half)
    local o=qInvRotate(qv,qw,origin-center)
    local d=qInvRotate(qv,qw,dir)
    local tmin,tmax=0,maxDist
    local hitAxis,hitSign=1,1
    local function slab(ov,dv,h,axis)
        if math.abs(dv)<EPS then
            return math.abs(ov)<=h
        end
        local inv=1/dv
        local t1=(-h-ov)*inv
        local t2=(h-ov)*inv
        local s=-1
        if t1>t2 then t1,t2=t2,t1;s=1 end
        if t1>tmin then tmin=t1;hitAxis=axis;hitSign=s end
        if t2<tmax then tmax=t2 end
        return tmin<=tmax
    end
    if not slab(o.X,d.X,half.X,1) or not slab(o.Y,d.Y,half.Y,2) or not slab(o.Z,d.Z,half.Z,3) then return nil end
    if tmin<0 or tmin>maxDist then return nil end
    local nl=if hitAxis==1 then vector.create(hitSign,0,0) elseif hitAxis==2 then vector.create(0,hitSign,0) else vector.create(0,0,hitSign)
    local p=origin+dir*tmin
    return tmin,qRotate(qv,qw,nl),p
end

function World:raycast(origin,direction,maxDistance,mask)
    maxDistance=maxDistance or math.huge
    mask=if mask==nil then 0xffffffff else mask
    local dl2=vector.dot(direction,direction)
    if dl2<=EPS then return nil end
    local dir=direction/math.sqrt(dl2)
    local bestT,bestI,bestN,bestP=maxDistance,nil,nil,nil
    for i=1,self.high do
        if not self.alive[i] or bit32.band(mask,self.group[i])==0 then continue end
        local s=self.shape[i]
        local t,n,p
        if s.kind==SHAPE_SPHERE then
            t,n,p=raySphere(origin,dir,bestT,self.pos[i],s.radius)
        elseif s.kind==SHAPE_BOX then
            t,n,p=rayBox(origin,dir,bestT,self.pos[i],self.qv[i],self.qw[i],s.halfExtents)
        else
            local a,b=capsuleEndpoints(self.pos[i],self.qv[i],self.qw[i],s.halfHeight)
            local m=(a+b)*0.5
            local ta,na,pa=raySphere(origin,dir,bestT,a,s.radius)
            local tm,nm,pm=raySphere(origin,dir,bestT,m,s.radius)
            local tb,nb,pb=raySphere(origin,dir,bestT,b,s.radius)
            t,n,p=ta,na,pa
            if tm and (not t or tm<t) then t,n,p=tm,nm,pm end
            if tb and (not t or tb<t) then t,n,p=tb,nb,pb end
        end
        if t and t<bestT then bestT,bestI,bestN,bestP=t,i,n,p end
    end
    if not bestI then return nil end
    return { handle=makeHandle(bestI,self.generation[bestI]), distance=bestT, point=bestP, normal=bestN, userData=self.userData[bestI] }
end

function World:contacts()
    local out=table.create(self.contactCount)
    for c=1,self.contactCount do
        local a,b=self.ca[c],self.cb[c]
        out[c]={
            a=makeHandle(a,self.generation[a]),
            b=makeHandle(b,self.generation[b]),
            normal=self.cn[c],point=self.cp[c],penetration=self.cpen[c],
            normalImpulse=self.cnormalImpulse[c],
        }
    end
    return out
end

function World:stats()
    local asleep=0
    for i=1,self.high do if self.alive[i] and self.sleeping[i] then asleep+=1 end end
    return {
        frame=self.frame,bodies=self.aliveCount,sleeping=asleep,
        candidates=self.candidateCount,contacts=self.contactCount,
        joints=(function() local n=0 for _ in self.joints do n+=1 end return n end)(),
    }
end

function M.newWorld(opts)
    opts=opts or {}
    local w=setmetatable({},World)
    w.gravity=opts.gravity or vector.create(0,-9.81,0)
    w.velocityIterations=opts.velocityIterations or 8
    w.positionCorrection=opts.positionCorrection or 0.65
    w.penetrationSlop=opts.penetrationSlop or 0.003
    w.restitutionThreshold=opts.restitutionThreshold or 1.0
    w.allowSleep=if opts.allowSleep==nil then true else opts.allowSleep
    w.sleepLinear=opts.sleepLinear or 0.04
    w.sleepAngular=opts.sleepAngular or 0.05
    w.sleepDelay=opts.sleepDelay or 0.6
    w.wakeSpeed=opts.wakeSpeed or 0.5
    w.wakePenetration=opts.wakePenetration or 0.05
    w.broadphase=opts.broadphase or "grid"
    assert(w.broadphase=="grid" or w.broadphase=="sap","physik.newWorld: broadphase must be grid or sap")
    w.broadphaseCellSize=opts.broadphaseCellSize or 1.0
    assert(w.broadphaseCellSize>0,"physik.newWorld: broadphaseCellSize must be > 0")
    w.gridHeads={};w.gridKeys={};w.gridNext={};w.gridLarge={};w.gridX={};w.gridY={};w.gridZ={};w.gridKey={}
    w.gridKeyX={};w.gridKeyY={};w.gridKeyZ={}
    w.frame=0
    w.high=0
    w.aliveCount=0
    w.free={}
    w.order={}
    w.generation={}
    w.alive={}
    w.bodyType={}
    w.shape={}
    w.pos={}
    w.qv={}
    w.qw={}
    w.axis0={};w.axis1={};w.axis2={}
    w.linVel={}
    w.angVel={}
    w.force={}
    w.torque={}
    w.invMass={}
    w.invInertia={}
    w.friction={}
    w.restitution={}
    w.linearDamping={}
    w.angularDamping={}
    w.group={}
    w.mask={}
    w.userData={}
    w.sleeping={}
    w.sleepTime={}
    w.aabbMin={}
    w.aabbMax={}
    w.contactCount=0
    w.candidateCount=0
    w.ca={};w.cb={};w.cn={};w.cp={};w.cpen={};w.cnormalImpulse={}
    w.ctangentImpulse1={};w.ctangentImpulse2={}
    w.cra={};w.crb={};w.cnormalMass={};w.cfriction={};w.crest={}
    w.ctangent1={};w.ctangent2={};w.ctangentMass1={};w.ctangentMass2={}
    w.joints={};w.nextJoint=1
    return w
end

M._constants = table.freeze({ sphere=SHAPE_SPHERE, box=SHAPE_BOX, capsule=SHAPE_CAPSULE })
return table.freeze(M)
