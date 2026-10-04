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
local JOINT_DISTANCE = 1
local JOINT_BALL = 2
local JOINT_HINGE = 3

-- Persistent contact cache entry layout. One array table per body pair, no child tables.
local CACHE_COUNT = 1
local CACHE_EPOCH = 2
local CACHE_GA = 3
local CACHE_GB = 4
local CACHE_MATCH_EPOCH = 5
local CACHE_USED_MASK = 6
local CACHE_HEADER = 6
local CACHE_STRIDE = 5 -- localA, localB, normal, normalImpulse, tangentImpulse

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

local function bodyPairKey(a,b)
    if b<a then a,b=b,a end
    return a*HANDLE_STRIDE+b
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

local function testBoxAxis(axis,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,feature,margin)
    local l2=vector.dot(axis,axis)
    if l2<=1e-10 then return true,bestPen,bestAxis,bestFeature end
    local n=if math.abs(l2-1)<1e-5 then axis else axis/math.sqrt(l2)
    local ra=ah.X*math.abs(vector.dot(n,a0))+ah.Y*math.abs(vector.dot(n,a1))+ah.Z*math.abs(vector.dot(n,a2))
    local rb=bh.X*math.abs(vector.dot(n,b0))+bh.Y*math.abs(vector.dot(n,b1))+bh.Z*math.abs(vector.dot(n,b2))
    local dn=vector.dot(d,n)
    local pen=ra+rb-math.abs(dn)
    if pen<=-(margin or 0) then return false,bestPen,bestAxis,bestFeature end
    if pen<bestPen then
        bestPen=pen
        bestAxis=if dn>=0 then n else -n
        bestFeature=feature
    end
    return true,bestPen,bestAxis,bestFeature
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

local function sphereSphere(pa, ra, pb, rb, margin)
    local d = pb - pa
    margin=margin or 0
    local rr = ra + rb
    local detect=rr+margin
    local d2 = vector.dot(d, d)
    if d2 >= detect*detect then return false end
    local dist = math.sqrt(math.max(d2, EPS))
    local n = if d2 > EPS then d / dist else YAXIS
    local penetration = rr - dist
    local point = pa + n * (ra - penetration * 0.5)
    return true, n, point, penetration
end

local function sphereBox(sp, radius, bp, bqv, bqw, half, margin)
    local rel = qInvRotate(bqv, bqw, sp - bp)
    local closest = vector.create(
        clamp(rel.X, -half.X, half.X),
        clamp(rel.Y, -half.Y, half.Y),
        clamp(rel.Z, -half.Z, half.Z)
    )
    local d = rel - closest -- box -> sphere in box-local space.
    local d2 = vector.dot(d, d)
    margin=margin or 0
    local detect=radius+margin
    if d2 > detect*detect then return false end
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

local function boxBox(ap,a0,a1,a2,ah,bp,b0,b1,b2,bh,margin)
    local d=bp-ap
    local bestPen,bestAxis,bestFeature=math.huge,nil,0
    local ok

    ok,bestPen,bestAxis,bestFeature=testBoxAxis(a0,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,1,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(a1,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,2,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(a2,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,3,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(b0,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,4,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(b1,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,5,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(b2,d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,6,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(vector.cross(a0,b0),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,7,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(vector.cross(a0,b1),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,8,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(vector.cross(a0,b2),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,9,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(vector.cross(a1,b0),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,10,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(vector.cross(a1,b1),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,11,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(vector.cross(a1,b2),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,12,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(vector.cross(a2,b0),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,13,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(vector.cross(a2,b1),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,14,margin); if not ok then return false end
    ok,bestPen,bestAxis,bestFeature=testBoxAxis(vector.cross(a2,b2),d,a0,a1,a2,ah,b0,b1,b2,bh,bestPen,bestAxis,bestFeature,15,margin); if not ok then return false end

    local n=bestAxis or safeUnit(d,YAXIS)
    local sa=supportBoxAxes(ap,a0,a1,a2,ah,n)
    local sb=supportBoxAxes(bp,b0,b1,b2,bh,-n)
    return true,n,(sa+sb)*0.5,bestPen,bestFeature
end

local function axisAt(a0,a1,a2,i)
    if i==1 then return a0 elseif i==2 then return a1 else return a2 end
end

local function extentAt(h,i)
    if i==1 then return h.X elseif i==2 then return h.Y else return h.Z end
end

local function otherAxes(a0,a1,a2,h,i)
    if i==1 then return a1,a2,h.Y,h.Z end
    if i==2 then return a0,a2,h.X,h.Z end
    return a0,a1,h.X,h.Y
end

local function clipPlane(input,count,output,origin,axis,limit,sign)
    if count<=0 then return 0 end
    local outCount=0
    local prev=input[count]
    local prevD=limit-sign*vector.dot(prev-origin,axis)
    local prevIn=prevD>=-1e-6
    for i=1,count do
        local cur=input[i]
        local curD=limit-sign*vector.dot(cur-origin,axis)
        local curIn=curD>=-1e-6
        if curIn~=prevIn then
            local denom=prevD-curD
            local t=if math.abs(denom)>EPS then prevD/denom else 0
            outCount+=1
            output[outCount]=prev+(cur-prev)*t
        end
        if curIn then
            outCount+=1
            output[outCount]=cur
        end
        prev,prevD,prevIn=cur,curD,curIn
    end
    return outCount
end

local function emitFaceManifold(self,a,b,n,feature,bestPen)
    local ref,inc,refIndex,refOut
    if feature<=3 then
        ref,inc,refIndex,refOut=a,b,feature,n
    else
        ref,inc,refIndex,refOut=b,a,feature-3,-n
    end

    local rp=self.pos[ref]
    local ra0,ra1,ra2=self.axis0[ref],self.axis1[ref],self.axis2[ref]
    local rh=self.shape[ref].halfExtents
    local refAxis=axisAt(ra0,ra1,ra2,refIndex)
    -- SAT oriented refOut already, but use its sign against the cached basis to land on the right face.
    local refSign=if vector.dot(refOut,refAxis)>=0 then 1 else -1
    refOut=refAxis*refSign
    local refCenter=rp+refOut*extentAt(rh,refIndex)
    local ru,rv,reu,rev=otherAxes(ra0,ra1,ra2,rh,refIndex)

    local ip=self.pos[inc]
    local ia0,ia1,ia2=self.axis0[inc],self.axis1[inc],self.axis2[inc]
    local ih=self.shape[inc].halfExtents
    local d0,d1,d2=vector.dot(refOut,ia0),vector.dot(refOut,ia1),vector.dot(refOut,ia2)
    local ad0,ad1,ad2=math.abs(d0),math.abs(d1),math.abs(d2)
    local incIndex=if ad0>=ad1 and ad0>=ad2 then 1 elseif ad1>=ad2 then 2 else 3
    local incAxis=axisAt(ia0,ia1,ia2,incIndex)
    local incDot=if incIndex==1 then d0 elseif incIndex==2 then d1 else d2
    local incOut=incAxis*(if incDot>0 then -1 else 1)
    local faceCenter=ip+incOut*extentAt(ih,incIndex)
    local iu,iv,ieu,iev=otherAxes(ia0,ia1,ia2,ih,incIndex)

    local ca,cb=self.clipA,self.clipB
    ca[1]=faceCenter+iu*ieu+iv*iev
    ca[2]=faceCenter-iu*ieu+iv*iev
    ca[3]=faceCenter-iu*ieu-iv*iev
    ca[4]=faceCenter+iu*ieu-iv*iev
    local count=4
    count=clipPlane(ca,count,cb,refCenter,ru,reu,1)
    count=clipPlane(cb,count,ca,refCenter,ru,reu,-1)
    count=clipPlane(ca,count,cb,refCenter,rv,rev,1)
    count=clipPlane(cb,count,ca,refCenter,rv,rev,-1)

    local points,depths=self.manifoldPoints,self.manifoldDepths
    local m=0
    for i=1,count do
        local p=ca[i]
        local sep=vector.dot(p-refCenter,refOut)
        if sep<=math.max(self.manifoldSlop,self.contactMargin) then
            m+=1
            local depth=math.max(-sep,0)
            points[m]=p-refOut*(sep*0.5)
            depths[m]=math.max(depth,math.min(bestPen,self.manifoldSlop))
        end
    end
    if m==0 then
        local sa=supportBoxAxes(self.pos[a],self.axis0[a],self.axis1[a],self.axis2[a],self.shape[a].halfExtents,n)
        local sb=supportBoxAxes(self.pos[b],self.axis0[b],self.axis1[b],self.axis2[b],self.shape[b].halfExtents,-n)
        self:_addContact(a,b,n,(sa+sb)*0.5,bestPen)
        return
    end

    if m<=4 then
        for i=1,m do self:_addContact(a,b,n,points[i],depths[i]) end
        return
    end

    -- Reduce an octagonal clipped face to four spatially spread points without allocating.
    local i1=1
    for i=2,m do if depths[i]>depths[i1] then i1=i end end
    local p1=points[i1]
    local i2=if i1==1 then 2 else 1
    local best=-1
    for i=1,m do
        if i~=i1 then
            local d=points[i]-p1
            local d2=vector.dot(d,d)
            if d2>best then best=d2;i2=i end
        end
    end
    local p2=points[i2]
    local i3=0;best=-1
    local edge=p2-p1
    for i=1,m do
        if i~=i1 and i~=i2 then
            local cr=vector.cross(edge,points[i]-p1)
            local area2=vector.dot(cr,cr)
            if area2>best then best=area2;i3=i end
        end
    end
    if i3==0 then i3=if i2~=1 and i1~=1 then 1 else 3 end
    local p3=points[i3]
    local i4=0;best=-1
    for i=1,m do
        if i~=i1 and i~=i2 and i~=i3 then
            local d1=points[i]-p1; local d2=points[i]-p2; local d3=points[i]-p3
            local score=math.min(vector.dot(d1,d1),vector.dot(d2,d2),vector.dot(d3,d3))
            if score>best then best=score;i4=i end
        end
    end
    self:_addContact(a,b,n,points[i1],depths[i1])
    self:_addContact(a,b,n,points[i2],depths[i2])
    self:_addContact(a,b,n,points[i3],depths[i3])
    if i4~=0 then self:_addContact(a,b,n,points[i4],depths[i4]) end
end

local function capsuleEndpoints(p, qv, qw, hh)
    local axis = qRotate(qv, qw, YAXIS)
    local d = axis * hh
    return p - d, p + d
end

local function capsuleSphere(cp, cqv, cqw, cr, chh, sp, sr, margin)
    local a,b = capsuleEndpoints(cp,cqv,cqw,chh)
    local q = closestPointSegment(sp,a,b)
    return sphereSphere(q,cr,sp,sr,margin)
end

local function capsuleCapsule(ap,aqv,aqw,ar,ahh,bp,bqv,bqw,br,bhh,margin)
    local a0,a1 = capsuleEndpoints(ap,aqv,aqw,ahh)
    local b0,b1 = capsuleEndpoints(bp,bqv,bqw,bhh)
    local ca,cb = closestSegmentSegment(a0,a1,b0,b1)
    return sphereSphere(ca,ar,cb,br,margin)
end

local function closestPointAabbLocal(p,h)
    return vector.create(clamp(p.X,-h.X,h.X),clamp(p.Y,-h.Y,h.Y),clamp(p.Z,-h.Z,h.Z))
end

local function segmentAabbInterval(a,b,h)
    local d=b-a
    local tmin,tmax=0,1
    local function slab(p,dp,lo,hi)
        if math.abs(dp)<=EPS then return p>=lo and p<=hi,tmin,tmax end
        local inv=1/dp
        local t1=(lo-p)*inv
        local t2=(hi-p)*inv
        if t1>t2 then t1,t2=t2,t1 end
        local nmin=math.max(tmin,t1)
        local nmax=math.min(tmax,t2)
        return nmin<=nmax,nmin,nmax
    end
    local ok
    ok,tmin,tmax=slab(a.X,d.X,-h.X,h.X); if not ok then return false end
    ok,tmin,tmax=slab(a.Y,d.Y,-h.Y,h.Y); if not ok then return false end
    ok,tmin,tmax=slab(a.Z,d.Z,-h.Z,h.Z); if not ok then return false end
    return true,tmin,tmax
end

local function considerSegmentEdge(a,b,e0,e1,bestD,bestS,bestB)
    local s,e=closestSegmentSegment(a,b,e0,e1)
    local dv=s-e
    local d2=vector.dot(dv,dv)
    if d2<bestD then return d2,s,e end
    return bestD,bestS,bestB
end

local function closestSegmentAabbLocal(a,b,h)
    local hit,t0,t1=segmentAabbInterval(a,b,h)
    if hit then
        local s=a+(b-a)*((t0+t1)*0.5)
        local dx=h.X-math.abs(s.X)
        local dy=h.Y-math.abs(s.Y)
        local dz=h.Z-math.abs(s.Z)
        if dx<=dy and dx<=dz then
            local out=vector.create(signnz(s.X),0,0)
            return s,s+out*dx,true,out,dx
        elseif dy<=dz then
            local out=vector.create(0,signnz(s.Y),0)
            return s,s+out*dy,true,out,dy
        else
            local out=vector.create(0,0,signnz(s.Z))
            return s,s+out*dz,true,out,dz
        end
    end

    local qa=closestPointAabbLocal(a,h)
    local da=a-qa
    local bestD=vector.dot(da,da)
    local bestS,bestB=a,qa
    local qb=closestPointAabbLocal(b,h)
    local db=b-qb
    local d2=vector.dot(db,db)
    if d2<bestD then bestD,bestS,bestB=d2,b,qb end

    local hx,hy,hz=h.X,h.Y,h.Z
    -- X-parallel edges.
    bestD,bestS,bestB=considerSegmentEdge(a,b,vector.create(-hx,-hy,-hz),vector.create(hx,-hy,-hz),bestD,bestS,bestB)
    bestD,bestS,bestB=considerSegmentEdge(a,b,vector.create(-hx,-hy, hz),vector.create(hx,-hy, hz),bestD,bestS,bestB)
    bestD,bestS,bestB=considerSegmentEdge(a,b,vector.create(-hx, hy,-hz),vector.create(hx, hy,-hz),bestD,bestS,bestB)
    bestD,bestS,bestB=considerSegmentEdge(a,b,vector.create(-hx, hy, hz),vector.create(hx, hy, hz),bestD,bestS,bestB)
    -- Y-parallel edges.
    bestD,bestS,bestB=considerSegmentEdge(a,b,vector.create(-hx,-hy,-hz),vector.create(-hx,hy,-hz),bestD,bestS,bestB)
    bestD,bestS,bestB=considerSegmentEdge(a,b,vector.create(-hx,-hy, hz),vector.create(-hx,hy, hz),bestD,bestS,bestB)
    bestD,bestS,bestB=considerSegmentEdge(a,b,vector.create( hx,-hy,-hz),vector.create( hx,hy,-hz),bestD,bestS,bestB)
    bestD,bestS,bestB=considerSegmentEdge(a,b,vector.create( hx,-hy, hz),vector.create( hx,hy, hz),bestD,bestS,bestB)
    -- Z-parallel edges.
    bestD,bestS,bestB=considerSegmentEdge(a,b,vector.create(-hx,-hy,-hz),vector.create(-hx,-hy,hz),bestD,bestS,bestB)
    bestD,bestS,bestB=considerSegmentEdge(a,b,vector.create(-hx, hy,-hz),vector.create(-hx, hy,hz),bestD,bestS,bestB)
    bestD,bestS,bestB=considerSegmentEdge(a,b,vector.create( hx,-hy,-hz),vector.create( hx,-hy,hz),bestD,bestS,bestB)
    bestD,bestS,bestB=considerSegmentEdge(a,b,vector.create( hx, hy,-hz),vector.create( hx, hy,hz),bestD,bestS,bestB)
    return bestS,bestB,false,nil,math.sqrt(bestD)
end

local function capsuleBox(cp,cqv,cqw,cr,chh,bp,bqv,bqw,bh,margin)
    local wa,wb=capsuleEndpoints(cp,cqv,cqw,chh)
    local a=qInvRotate(bqv,bqw,wa-bp)
    local b=qInvRotate(bqv,bqw,wb-bp)
    local spine,boxPoint,inside,outward,measure=closestSegmentAabbLocal(a,b,bh)
    margin=margin or 0
    local nLocal,penetration
    if inside then
        nLocal=-outward
        penetration=cr+measure
    else
        local delta=spine-boxPoint -- box -> capsule spine
        local d2=vector.dot(delta,delta)
        local detect=cr+margin
        if d2>detect*detect then return false end
        local dist=math.sqrt(math.max(d2,0))
        if dist>EPS then nLocal=-(delta/dist) else nLocal=YAXIS end
        penetration=cr-dist
    end
    local normal=qRotate(bqv,bqw,nLocal)
    local point=bp+qRotate(bqv,bqw,boxPoint)
    return true,normal,point,penetration
end

local function sweepSphereSphere(a,b,radius,center,targetRadius)
    local d=b-a
    local m=a-center
    local rr=radius+targetRadius
    local c=vector.dot(m,m)-rr*rr
    if c<=0 then return nil end -- already overlapping: leave this to the discrete solver.
    local aa=vector.dot(d,d)
    if aa<=EPS then return nil end
    local bb=vector.dot(m,d)
    if bb>=0 then return nil end
    local disc=bb*bb-aa*c
    if disc<0 then return nil end
    local t=(-bb-math.sqrt(disc))/aa
    if t<0 or t>1 then return nil end
    local hit=a+d*t
    return t,safeUnit(hit-center,YAXIS)
end

local function sweepSphereBox(a,b,radius,bp,bqv,bqw,bh)
    local la=qInvRotate(bqv,bqw,a-bp)
    local lb=qInvRotate(bqv,bqw,b-bp)
    local h=vector.create(bh.X+radius,bh.Y+radius,bh.Z+radius)
    local inside=math.abs(la.X)<=h.X and math.abs(la.Y)<=h.Y and math.abs(la.Z)<=h.Z
    if inside then return nil end
    local hit,t0=segmentAabbInterval(la,lb,h)
    if not hit or t0<0 or t0>1 then return nil end
    local p=la+(lb-la)*t0
    local dx=math.abs(math.abs(p.X)-h.X)
    local dy=math.abs(math.abs(p.Y)-h.Y)
    local dz=math.abs(math.abs(p.Z)-h.Z)
    local nLocal
    if dx<=dy and dx<=dz then nLocal=vector.create(signnz(p.X),0,0)
    elseif dy<=dz then nLocal=vector.create(0,signnz(p.Y),0)
    else nLocal=vector.create(0,0,signnz(p.Z)) end
    return t0,qRotate(bqv,bqw,nLocal)
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
    self.ccd[i] = desc.ccd == true
    self.prevPos[i] = self.pos[i]
    self.sleeping[i] = false
    self.sleepTime[i] = 0
    self.aabbMin[i] = ZERO
    self.aabbMax[i] = ZERO
    return makeHandle(i,self.generation[i])
end

function World:removeBody(handle)
    local i = self:_check(handle)
    for id,j in self.joints do
        if j.a==i or j.b==i then self:_dropJoint(id,j) end
    end
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
        kind=JOINT_DISTANCE,
        a=a,b=b, ga=self.generation[a],gb=self.generation[b],
        la=la,lb=lb,
        rest=opts.length or vector.magnitude(wb-wa),
        stiffness=opts.stiffness or 0.25,
        damping=opts.damping or 0.05,
        impulse=0,
        collideConnected=opts.collideConnected==true,
    }
    if not self.joints[id].collideConnected then
        local key=bodyPairKey(a,b)
        self.jointNoCollide[key]=(self.jointNoCollide[key] or 0)+1
    end
    return id
end

function World:addBallJoint(aHandle,bHandle,opts)
    local a,b=self:_check(aHandle),self:_check(bHandle)
    opts=opts or {}
    local id=self.nextJoint
    self.nextJoint+=1
    self.joints[id]={kind=JOINT_BALL,a=a,b=b,ga=self.generation[a],gb=self.generation[b],la=opts.localAnchorA or ZERO,lb=opts.localAnchorB or ZERO,stiffness=opts.stiffness or 0.25,damping=opts.damping or 0.05,collideConnected=opts.collideConnected==true}
    if not self.joints[id].collideConnected then
        local key=bodyPairKey(a,b)
        self.jointNoCollide[key]=(self.jointNoCollide[key] or 0)+1
    end
    return id
end

function World:addHingeJoint(aHandle,bHandle,opts)
    local a,b=self:_check(aHandle),self:_check(bHandle)
    opts=opts or {}
    local id=self.nextJoint
    self.nextJoint+=1
    self.joints[id]={kind=JOINT_HINGE,a=a,b=b,ga=self.generation[a],gb=self.generation[b],la=opts.localAnchorA or ZERO,lb=opts.localAnchorB or ZERO,axisA=safeUnit(opts.localAxisA or YAXIS,YAXIS),axisB=safeUnit(opts.localAxisB or YAXIS,YAXIS),stiffness=opts.stiffness or 0.25,damping=opts.damping or 0.05,angularStiffness=opts.angularStiffness or 0.2,motorSpeed=opts.motorSpeed,maxMotorTorque=opts.maxMotorTorque or math.huge,motorImpulse=0,collideConnected=opts.collideConnected==true}
    if not self.joints[id].collideConnected then
        local key=bodyPairKey(a,b)
        self.jointNoCollide[key]=(self.jointNoCollide[key] or 0)+1
    end
    return id
end

function World:_dropJoint(id,j)
    j=j or self.joints[id]
    if j==nil then return end
    if not j.collideConnected then
        local key=bodyPairKey(j.a,j.b)
        local count=self.jointNoCollide[key] or 0
        if count<=1 then self.jointNoCollide[key]=nil else self.jointNoCollide[key]=count-1 end
    end
    self.joints[id]=nil
end

function World:removeJoint(id)
    self:_dropJoint(id)
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
    local margin=self.contactMargin
    if margin>0 then
        local skin=vector.create(margin,margin,margin)
        self.aabbMin[i]-=skin
        self.aabbMax[i]+=skin
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
local applyPairImpulse

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
        ok,n,p,pen=sphereSphere(ap,sa.radius,bp,sb.radius,self.contactMargin)
    elseif sa.kind==SHAPE_SPHERE and sb.kind==SHAPE_BOX then
        ok,n,p,pen=sphereBox(ap,sa.radius,bp,bqv,bqw,sb.halfExtents,self.contactMargin)
    elseif sa.kind==SHAPE_BOX and sb.kind==SHAPE_SPHERE then
        ok,n,p,pen=sphereBox(bp,sb.radius,ap,aqv,aqw,sa.halfExtents,self.contactMargin)
        if ok then n=-n end
    elseif sa.kind==SHAPE_BOX and sb.kind==SHAPE_BOX then
        local feature
        ok,n,p,pen,feature=boxBox(ap,self.axis0[a],self.axis1[a],self.axis2[a],sa.halfExtents,bp,self.axis0[b],self.axis1[b],self.axis2[b],sb.halfExtents,self.contactMargin)
        if ok and feature<=6 then
            emitFaceManifold(self,a,b,n,feature,pen)
            return
        end
    elseif sa.kind==SHAPE_CAPSULE and sb.kind==SHAPE_SPHERE then
        ok,n,p,pen=capsuleSphere(ap,aqv,aqw,sa.radius,sa.halfHeight,bp,sb.radius,self.contactMargin)
    elseif sa.kind==SHAPE_SPHERE and sb.kind==SHAPE_CAPSULE then
        ok,n,p,pen=capsuleSphere(bp,bqv,bqw,sb.radius,sb.halfHeight,ap,sa.radius,self.contactMargin)
        if ok then n=-n end
    elseif sa.kind==SHAPE_CAPSULE and sb.kind==SHAPE_CAPSULE then
        ok,n,p,pen=capsuleCapsule(ap,aqv,aqw,sa.radius,sa.halfHeight,bp,bqv,bqw,sb.radius,sb.halfHeight,self.contactMargin)
    elseif sa.kind==SHAPE_CAPSULE and sb.kind==SHAPE_BOX then
        ok,n,p,pen=capsuleBox(ap,aqv,aqw,sa.radius,sa.halfHeight,bp,bqv,bqw,sb.halfExtents,self.contactMargin)
    elseif sa.kind==SHAPE_BOX and sb.kind==SHAPE_CAPSULE then
        ok,n,p,pen=capsuleBox(bp,bqv,bqw,sb.radius,sb.halfHeight,ap,aqv,aqw,sa.halfExtents,self.contactMargin)
        if ok then n=-n end
    end
    if ok then self:_addContact(a,b,n,p,pen) end
end

function World:_tryPair(a,b)
    if a==b or not self.alive[a] or not self.alive[b] then return end
    if b<a then a,b=b,a end -- stable pair orientation for persistent contact matching.
    if self.jointNoCollide[bodyPairKey(a,b)] then return end
    if self.bodyType[a]==BODY_STATIC and self.bodyType[b]==BODY_STATIC then return end
    if bit32.band(self.mask[a],self.group[b])==0 or bit32.band(self.mask[b],self.group[a])==0 then return end
    local amin,amax=self.aabbMin[a],self.aabbMax[a]
    local bmin,bmax=self.aabbMin[b],self.aabbMax[b]
    if bmin.X>amax.X or bmax.X<amin.X or bmin.Y>amax.Y or bmax.Y<amin.Y or bmin.Z>amax.Z or bmax.Z<amin.Z then return end
    self.candidateCount+=1
    local before=self.contactCount
    self:_narrow(a,b)
    if self.contactCount>before then self.manifoldCount+=1 end
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
        if #large<=self.directLargeThreshold then
            -- A handful of oversized bodies (the usual case is one world/floor collider)
            -- should not force an all-body sweep that merely rediscovers small-small pairs.
            -- Emit each large-small pair once, plus large-large pairs by stable body index.
            for li=1,#large do
                local a=large[li]
                for b=1,self.high do
                    if self.alive[b] and b~=a and (self.gridKey[b]~=nil or b>a) then
                        self:_tryPair(a,b)
                    end
                end
            end
        else
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
end

function World:_collide()
    self.contactCount=0
    self.candidateCount=0
    self.manifoldCount=0
    self.cacheEpoch+=1
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

local function contactAxisResponse(self,i,r,axis)
    if self.invMass[i]==0 or self.sleeping[i] then return 0,ZERO end
    local rx=vector.cross(r,axis)
    local angular=if vector.dot(rx,rx)>EPS then applyInvInertia(self,i,rx) else ZERO
    return self.invMass[i]+vector.dot(axis,vector.cross(angular,r)),angular
end

local function pairKey(a,b)
    return a*HANDLE_STRIDE+b
end

local function unionRoot(parent,x)
    local p=parent[x]
    while p~=x do
        local gp=parent[p]
        parent[x]=gp
        x=p
        p=gp
    end
    return x
end

local function unionBodies(parent,size,a,b)
    local ra,rb=unionRoot(parent,a),unionRoot(parent,b)
    if ra==rb then return ra end
    if size[ra]<size[rb] then ra,rb=rb,ra end
    parent[rb]=ra
    size[ra]+=size[rb]
    return ra
end

function World:_buildIslands()
    local parent,size=self.islandParent,self.islandSize
    for i=1,self.high do
        if self.alive[i] and self.bodyType[i]==BODY_DYNAMIC then
            parent[i]=i
            size[i]=1
        else
            parent[i]=0
            size[i]=0
        end
    end
    for c=1,self.contactCount do
        local a,b=self.ca[c],self.cb[c]
        if self.bodyType[a]==BODY_DYNAMIC and self.bodyType[b]==BODY_DYNAMIC then
            unionBodies(parent,size,a,b)
        end
    end
    for _,j in self.joints do
        local a,b=j.a,j.b
        if self.alive[a] and self.alive[b] and self.generation[a]==j.ga and self.generation[b]==j.gb
            and self.bodyType[a]==BODY_DYNAMIC and self.bodyType[b]==BODY_DYNAMIC then
            unionBodies(parent,size,a,b)
        end
    end

    local awake=self.islandAwake
    local roots=self.islandRoots
    for k=#roots,1,-1 do roots[k]=nil end
    for i=1,self.high do awake[i]=false end
    local count=0
    for i=1,self.high do
        if parent[i]~=0 then
            local r=unionRoot(parent,i)
            parent[i]=r
            if r==i then count+=1; roots[count]=r end
            if not self.sleeping[i] then awake[r]=true end
        end
    end
    -- Wake is an island property: an explicit impulse or meaningful impact on one member
    -- wakes every dynamically-connected member before solving.
    for i=1,self.high do
        local r=parent[i]
        if r~=0 and awake[r] then
            self.sleeping[i]=false
        end
    end
    self.islandCount=count
end

function World:_matchCachedContact(c,t1,t2)
    if not self.warmStarting or self.cbias[c] > 0 then return end
    local a,b=self.ca[c],self.cb[c]
    local e=self.contactCache[pairKey(a,b)]
    if e==nil or e[CACHE_EPOCH]~=self.cacheEpoch-1 or e[CACHE_GA]~=self.generation[a] or e[CACHE_GB]~=self.generation[b] then return end
    local p=self.cp[c]
    local la=qInvRotate(self.qv[a],self.qw[a],p-self.pos[a])
    local lb=qInvRotate(self.qv[b],self.qw[b],p-self.pos[b])
    self.cla[c],self.clb[c]=la,lb
    if e[CACHE_MATCH_EPOCH]~=self.cacheEpoch then
        e[CACHE_MATCH_EPOCH]=self.cacheEpoch
        e[CACHE_USED_MASK]=0
    end
    local best,bestD=0,self.manifoldMatchDistance*self.manifoldMatchDistance*2
    local n=self.cn[c]
    local used=e[CACHE_USED_MASK] or 0
    local count=e[CACHE_COUNT] or 0
    for j=1,count do
        local bit=bit32.lshift(1,j-1)
        local k=CACHE_HEADER+(j-1)*CACHE_STRIDE
        if bit32.band(used,bit)==0 and vector.dot(e[k+3],n)>0.85 then
            local da=la-e[k+1]
            local db=lb-e[k+2]
            local d=vector.dot(da,da)+vector.dot(db,db)
            if d<bestD then bestD=d;best=j end
        end
    end
    if best~=0 then
        local bit=bit32.lshift(1,best-1)
        e[CACHE_USED_MASK]=bit32.bor(used,bit)
        local k=CACHE_HEADER+(best-1)*CACHE_STRIDE
        local factor=self.warmStartFactor
        self.cnormalImpulse[c]=e[k+4]*factor
        local tangent=e[k+5]*factor
        self.ctangentImpulse1[c]=vector.dot(tangent,t1)
        self.ctangentImpulse2[c]=vector.dot(tangent,t2)
        self.warmStartedCount+=1
    end
end

function World:_storeContactCache()
    local epoch=self.cacheEpoch
    local active=self.activeContacts
    for ci=1,self.activeContactCount do
        local c=active[ci]
        local a,b=self.ca[c],self.cb[c]
        local key=pairKey(a,b)
        local e=self.contactCache[key]
        if e==nil then
            e=table.create(CACHE_HEADER+4*CACHE_STRIDE)
            e[CACHE_COUNT]=0;e[CACHE_EPOCH]=0;e[CACHE_GA]=0;e[CACHE_GB]=0;e[CACHE_MATCH_EPOCH]=0;e[CACHE_USED_MASK]=0
            self.contactCache[key]=e
        end
        if e[CACHE_EPOCH]~=epoch then
            e[CACHE_EPOCH]=epoch
            e[CACHE_COUNT]=0
            e[CACHE_GA]=self.generation[a]
            e[CACHE_GB]=self.generation[b]
        end
        local count=e[CACHE_COUNT]
        if count<4 then
            local j=count+1
            e[CACHE_COUNT]=j
            local k=CACHE_HEADER+(j-1)*CACHE_STRIDE
            local p=self.cp[c]
            e[k+1]=self.cla[c] or qInvRotate(self.qv[a],self.qw[a],p-self.pos[a])
            e[k+2]=self.clb[c] or qInvRotate(self.qv[b],self.qw[b],p-self.pos[b])
            e[k+3]=self.cn[c]
            e[k+4]=self.cnormalImpulse[c]
            e[k+5]=self.ctangent1[c]*self.ctangentImpulse1[c]+self.ctangent2[c]*self.ctangentImpulse2[c]
        end
    end
    if epoch%64==0 then
        for key,e in self.contactCache do
            if epoch-e[CACHE_EPOCH]>2 then self.contactCache[key]=nil end
        end
    end
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

function World:_prepareContacts(dt)
    self.warmStartedCount=0
    local activeList=self.activeContacts
    local activeCount=0
    for c=1,self.contactCount do
        local a,b=self.ca[c],self.cb[c]
        local active=(self.invMass[a]>0 and not self.sleeping[a]) or (self.invMass[b]>0 and not self.sleeping[b])
        self.cactive[c]=active
        if active then
            activeCount+=1
            activeList[activeCount]=c
        end
        if not active then
            self.cnormalMass[c]=0
            self.ctangentMass1[c]=0
            self.ctangentMass2[c]=0
            continue
        end
        local p,n=self.cp[c],self.cn[c]
        local ra,rb=p-self.pos[a],p-self.pos[b]
        self.cra[c],self.crb[c]=ra,rb
        self.cla[c]=qInvRotate(self.qv[a],self.qw[a],ra)
        self.clb[c]=qInvRotate(self.qv[b],self.qw[b],rb)
        local da,naA=contactAxisResponse(self,a,ra,n)
        local db,naB=contactAxisResponse(self,b,rb,n)
        local d=da+db
        self.cnormalMass[c]=if d>EPS then 1/d else 0
        self.cnormalAngA[c],self.cnormalAngB[c]=naA,naB
        local t1,t2=contactTangents(n)
        self.ctangent1[c],self.ctangent2[c]=t1,t2
        local d1a,t1a=contactAxisResponse(self,a,ra,t1)
        local d1b,t1b=contactAxisResponse(self,b,rb,t1)
        local d2a,t2a=contactAxisResponse(self,a,ra,t2)
        local d2b,t2b=contactAxisResponse(self,b,rb,t2)
        self.ctangentMass1[c]=if d1a+d1b>EPS then 1/(d1a+d1b) else 0
        self.ctangentMass2[c]=if d2a+d2b>EPS then 1/(d2a+d2b) else 0
        self.ctangentAng1A[c],self.ctangentAng1B[c]=t1a,t1b
        self.ctangentAng2A[c],self.ctangentAng2B[c]=t2a,t2b
        self.cfriction[c]=math.sqrt(self.friction[a]*self.friction[b])
        self.crest[c]=math.max(self.restitution[a],self.restitution[b])
        local biasSlop=math.max(self.penetrationSlop,self.baumgarteSlop)
        local depth=math.max(self.cpen[c]-biasSlop,0)
        self.cbias[c]=math.min(depth*self.baumgarte/math.max(dt,1e-5),self.maxDepenetrationVelocity)
        self:_matchCachedContact(c,t1,t2)
    end
    for i=#activeList,activeCount+1,-1 do activeList[i]=nil end
    self.activeContactCount=activeCount
    -- Seed this frame's solve with last frame's accumulated impulses. This is done after
    -- every contact mass is prepared, so a warm start never changes the effective mass basis.
    for ci=1,activeCount do
        local c=activeList[ci]
        local ni=self.cnormalImpulse[c]
        local ti1,ti2=self.ctangentImpulse1[c],self.ctangentImpulse2[c]
        if ni~=0 or ti1~=0 or ti2~=0 then
            local a,b=self.ca[c],self.cb[c]
            local linear=self.cn[c]*ni+self.ctangent1[c]*ti1+self.ctangent2[c]*ti2
            if self.invMass[a]>0 and not self.sleeping[a] then
                self.linVel[a]-=linear*self.invMass[a]
                self.angVel[a]-=self.cnormalAngA[c]*ni+self.ctangentAng1A[c]*ti1+self.ctangentAng2A[c]*ti2
            end
            if self.invMass[b]>0 and not self.sleeping[b] then
                self.linVel[b]+=linear*self.invMass[b]
                self.angVel[b]+=self.cnormalAngB[c]*ni+self.ctangentAng1B[c]*ti1+self.ctangentAng2B[c]*ti2
            end
        end
    end
end

applyPairImpulse = function(self,a,b,ra,rb,impulse)
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
    local active=self.activeContacts
    for iter=1,iterations do
        for ci=1,self.activeContactCount do
            local c=active[ci]
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
                local target=math.max(bounce,self.cbias[c])
                local lambda=(target-vn)*normalMass
                local old=self.cnormalImpulse[c]
                local now=math.max(old+lambda,0)
                lambda=now-old
                self.cnormalImpulse[c]=now
                if lambda~=0 then
                    if self.invMass[a]>0 and not self.sleeping[a] then
                        self.linVel[a]-=n*(self.invMass[a]*lambda)
                        self.angVel[a]-=self.cnormalAngA[c]*lambda
                    end
                    if self.invMass[b]>0 and not self.sleeping[b] then
                        self.linVel[b]+=n*(self.invMass[b]*lambda)
                        self.angVel[b]+=self.cnormalAngB[c]*lambda
                    end
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
                if d1~=0 or d2~=0 then
                    local linear=t1*d1+t2*d2
                    if self.invMass[a]>0 and not self.sleeping[a] then
                        self.linVel[a]-=linear*self.invMass[a]
                        self.angVel[a]-=self.ctangentAng1A[c]*d1+self.ctangentAng2A[c]*d2
                    end
                    if self.invMass[b]>0 and not self.sleeping[b] then
                        self.linVel[b]+=linear*self.invMass[b]
                        self.angVel[b]+=self.ctangentAng1B[c]*d1+self.ctangentAng2B[c]*d2
                    end
                end
            end
        end
        self:_solveJoints(dt)
    end
end

local function solvePointAxis(self,a,b,ra,rb,error,axis,stiffness,damping,dt)
    local denom=impulseDenom(self,a,ra,axis)+impulseDenom(self,b,rb,axis)
    if denom<=EPS then return end
    local rv=pointVelocity(self,b,rb)-pointVelocity(self,a,ra)
    local cdot=vector.dot(rv,axis)
    local bias=vector.dot(error,axis)*stiffness/math.max(dt,1e-5)
    local lambda=-(cdot+bias+damping*cdot)/denom
    applyPairImpulse(self,a,b,ra,rb,axis*lambda)
end

local function angularDenom(self,i,axis)
    if self.invMass[i]==0 or self.sleeping[i] then return 0 end
    return vector.dot(axis,applyInvInertia(self,i,axis))
end

local function applyAngularPairImpulse(self,a,b,impulse)
    if self.invMass[a]>0 and not self.sleeping[a] then self.angVel[a]-=applyInvInertia(self,a,impulse) end
    if self.invMass[b]>0 and not self.sleeping[b] then self.angVel[b]+=applyInvInertia(self,b,impulse) end
end

local function solveBallAxes(self,j,a,b,ra,rb,wa,wb,dt)
    local error=wb-wa
    solvePointAxis(self,a,b,ra,rb,error,XAXIS,j.stiffness,j.damping,dt)
    solvePointAxis(self,a,b,ra,rb,error,YAXIS,j.stiffness,j.damping,dt)
    solvePointAxis(self,a,b,ra,rb,error,ZAXIS,j.stiffness,j.damping,dt)
end

function World:_solveJoints(dt)
    for id,j in self.joints do
        local a,b=j.a,j.b
        if not self.alive[a] or not self.alive[b] or self.generation[a]~=j.ga or self.generation[b]~=j.gb then self:_dropJoint(id,j); continue end
        if (self.invMass[a]==0 or self.sleeping[a]) and (self.invMass[b]==0 or self.sleeping[b]) then continue end
        local ra=qRotate(self.qv[a],self.qw[a],j.la)
        local rb=qRotate(self.qv[b],self.qw[b],j.lb)
        local wa=self.pos[a]+ra
        local wb=self.pos[b]+rb
        if j.kind==JOINT_DISTANCE then
            local d=wb-wa
            local len=vector.magnitude(d)
            local n=if len>EPS then d/len else YAXIS
            local rv=pointVelocity(self,b,rb)-pointVelocity(self,a,ra)
            local denom=impulseDenom(self,a,ra,n)+impulseDenom(self,b,rb,n)
            if denom>EPS then
                local error=len-j.rest
                local cdot=vector.dot(rv,n)
                local bias=(error*j.stiffness)/math.max(dt,1e-5)
                local lambda=-(cdot+bias+j.damping*cdot)/denom
                applyPairImpulse(self,a,b,ra,rb,n*lambda)
            end
        else
            solveBallAxes(self,j,a,b,ra,rb,wa,wb,dt)
            if j.kind==JOINT_HINGE then
                local axisA=qRotate(self.qv[a],self.qw[a],j.axisA)
                local axisB=qRotate(self.qv[b],self.qw[b],j.axisB)
                local hinge=safeUnit(axisA+axisB,axisA)
                local t1,t2=contactTangents(hinge)
                local error=vector.cross(axisA,axisB)
                local wrel=self.angVel[b]-self.angVel[a]
                local stiffness=j.angularStiffness/math.max(dt,1e-5)
                local denom1=angularDenom(self,a,t1)+angularDenom(self,b,t1)
                if denom1>EPS then
                    local lambda=-(vector.dot(wrel,t1)+vector.dot(error,t1)*stiffness)/denom1
                    applyAngularPairImpulse(self,a,b,t1*lambda)
                end
                wrel=self.angVel[b]-self.angVel[a]
                local denom2=angularDenom(self,a,t2)+angularDenom(self,b,t2)
                if denom2>EPS then
                    local lambda=-(vector.dot(wrel,t2)+vector.dot(error,t2)*stiffness)/denom2
                    applyAngularPairImpulse(self,a,b,t2*lambda)
                end
                if j.motorSpeed~=nil and j.maxMotorTorque>0 then
                    wrel=self.angVel[b]-self.angVel[a]
                    local denom=angularDenom(self,a,hinge)+angularDenom(self,b,hinge)
                    if denom>EPS then
                        local lambda=(j.motorSpeed-vector.dot(wrel,hinge))/denom
                        local old=j.motorImpulse
                        local limit=j.maxMotorTorque*dt
                        local now=clamp(old+lambda,-limit,limit)
                        lambda=now-old
                        j.motorImpulse=now
                        applyAngularPairImpulse(self,a,b,hinge*lambda)
                    end
                end
            end
        end
    end
end

function World:_correctPositions()
    local slop=self.penetrationSlop
    local iterations=self.positionIterations
    local percent=self.positionCorrection/iterations
    local active=self.activeContacts
    for _=1,iterations do
        for ci=1,self.activeContactCount do
            local c=active[ci]
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
            if self.ccd[i] then self.prevPos[i]=self.pos[i] end
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

function World:_ccd()
    self.ccdHitCount=0
    for i=1,self.high do
        if not self.alive[i] or self.bodyType[i]~=BODY_DYNAMIC or self.sleeping[i] or not self.ccd[i] then continue end
        local shape=self.shape[i]
        if shape.kind~=SHAPE_SPHERE then continue end -- v0.2 CCD is intentionally sphere-first.
        local start=self.prevPos[i]
        local finish=self.pos[i]
        local motion=finish-start
        local len2=vector.dot(motion,motion)
        if len2<=EPS then continue end
        local bestT,bestN,bestJ=1,nil,nil
        for j=1,self.high do
            if j==i or not self.alive[j] or self.bodyType[j]==BODY_DYNAMIC then continue end
            if bit32.band(self.mask[i],self.group[j])==0 or bit32.band(self.mask[j],self.group[i])==0 then continue end
            local target=self.shape[j]
            local t,n
            if target.kind==SHAPE_SPHERE then
                t,n=sweepSphereSphere(start,finish,shape.radius,self.pos[j],target.radius)
            elseif target.kind==SHAPE_BOX then
                t,n=sweepSphereBox(start,finish,shape.radius,self.pos[j],self.qv[j],self.qw[j],target.halfExtents)
            end
            if t and t<bestT then bestT,bestN,bestJ=t,n,j end
        end
        if bestN then
            local len=math.sqrt(len2)
            local backoff=if len>EPS then math.min(self.ccdSlop/len,bestT) else 0
            self.pos[i]=start+motion*(bestT-backoff)
            local vn=vector.dot(self.linVel[i],bestN)
            if vn<0 then
                local rest=math.max(self.restitution[i],self.restitution[bestJ])
                self.linVel[i]-=bestN*(vn*(1+rest))
            end
            self.ccdHitCount+=1
        end
    end
end

function World:_sleep(dt)
    if not self.allowSleep then return end
    local parent=self.islandParent
    local canSleep,minTime=self.islandCanSleep,self.islandMinTime
    local lin2=self.sleepLinear*self.sleepLinear
    local ang2=self.sleepAngular*self.sleepAngular
    for k=1,#self.islandRoots do
        local r=self.islandRoots[k]
        canSleep[r]=true
        minTime[r]=math.huge
    end
    for i=1,self.high do
        local r=parent[i]
        if r~=0 then
            if vector.dot(self.linVel[i],self.linVel[i])>=lin2 or vector.dot(self.angVel[i],self.angVel[i])>=ang2 then
                canSleep[r]=false
            end
        end
    end
    for i=1,self.high do
        local r=parent[i]
        if r~=0 then
            if canSleep[r] then
                self.sleepTime[i]+=dt
                if self.sleepTime[i]<minTime[r] then minTime[r]=self.sleepTime[i] end
            else
                self.sleepTime[i]=0
                self.sleeping[i]=false
            end
        end
    end
    for i=1,self.high do
        local r=parent[i]
        if r~=0 and canSleep[r] and minTime[r]>=self.sleepDelay then
            self.sleeping[i]=true
            self.linVel[i],self.angVel[i]=ZERO,ZERO
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
        self:_ccd()
        self:_collide()
        self:_wakeContacts()
        self:_buildIslands()
        self:_prepareContacts(h)
        self:_solveContacts(h)
        self:_storeContactCache()
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
        frame=self.frame,bodies=self.aliveCount,sleeping=asleep,islands=self.islandCount,
        candidates=self.candidateCount,contacts=self.contactCount,activeContacts=self.activeContactCount,manifolds=self.manifoldCount,warmStarted=self.warmStartedCount,ccdHits=self.ccdHitCount,
        joints=(function() local n=0 for _ in self.joints do n+=1 end return n end)(),
    }
end

function M.newWorld(opts)
    opts=opts or {}
    local w=setmetatable({},World)
    w.gravity=opts.gravity or vector.create(0,-9.81,0)
    w.velocityIterations=opts.velocityIterations or 4
    w.positionIterations=opts.positionIterations or 1
    assert(w.positionIterations>=1 and w.positionIterations%1==0,"physik.newWorld: positionIterations must be a positive integer")
    w.positionCorrection=if opts.positionCorrection==nil then 0.12 else opts.positionCorrection
    w.baumgarte=if opts.baumgarte==nil then 0.03 else opts.baumgarte
    w.baumgarteSlop=if opts.baumgarteSlop==nil then 0.05 else opts.baumgarteSlop
    w.maxDepenetrationVelocity=opts.maxDepenetrationVelocity or 1.0
    w.penetrationSlop=opts.penetrationSlop or 0.003
    w.restitutionThreshold=opts.restitutionThreshold or 1.0
    w.allowSleep=if opts.allowSleep==nil then true else opts.allowSleep
    w.sleepLinear=opts.sleepLinear or 0.04
    w.sleepAngular=opts.sleepAngular or 0.05
    w.sleepDelay=opts.sleepDelay or 0.6
    w.wakeSpeed=opts.wakeSpeed or 0.5
    w.wakePenetration=opts.wakePenetration or 0.05
    w.warmStarting=if opts.warmStarting==nil then true else opts.warmStarting
    w.warmStartFactor=opts.warmStartFactor or 0.85
    w.manifoldMatchDistance=opts.manifoldMatchDistance or 0.08
    w.contactMargin=if opts.contactMargin==nil then 0.015 else opts.contactMargin
    w.manifoldSlop=opts.manifoldSlop or w.contactMargin
    w.ccdSlop=opts.ccdSlop or 0.001
    w.broadphase=opts.broadphase or "grid"
    assert(w.broadphase=="grid" or w.broadphase=="sap","physik.newWorld: broadphase must be grid or sap")
    w.broadphaseCellSize=opts.broadphaseCellSize or 1.0
    assert(w.broadphaseCellSize>0,"physik.newWorld: broadphaseCellSize must be > 0")
    w.directLargeThreshold=opts.directLargeThreshold or 16
    assert(w.directLargeThreshold>=0 and w.directLargeThreshold%1==0,"physik.newWorld: directLargeThreshold must be a non-negative integer")
    w.gridHeads={};w.gridKeys={};w.gridNext={};w.gridLarge={};w.gridX={};w.gridY={};w.gridZ={};w.gridKey={}
    w.gridKeyX={};w.gridKeyY={};w.gridKeyZ={}
    w.frame=0
    w.cacheEpoch=0
    w.ccdHitCount=0
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
    w.ccd={};w.prevPos={}
    w.sleeping={}
    w.sleepTime={}
    w.aabbMin={}
    w.aabbMax={}
    w.contactCount=0
    w.candidateCount=0
    w.manifoldCount=0
    w.warmStartedCount=0
    w.activeContactCount=0
    w.islandCount=0
    w.ca={};w.cb={};w.cn={};w.cp={};w.cpen={};w.cactive={};w.cnormalImpulse={}
    w.ctangentImpulse1={};w.ctangentImpulse2={}
    w.cra={};w.crb={};w.cla={};w.clb={};w.cnormalMass={};w.cfriction={};w.crest={};w.cbias={}
    w.cnormalAngA={};w.cnormalAngB={}
    w.ctangent1={};w.ctangent2={};w.ctangentMass1={};w.ctangentMass2={}
    w.ctangentAng1A={};w.ctangentAng1B={};w.ctangentAng2A={};w.ctangentAng2B={}
    w.contactCache={}
    w.activeContacts={}
    w.clipA={};w.clipB={};w.manifoldPoints={};w.manifoldDepths={}
    w.islandParent={};w.islandSize={};w.islandAwake={};w.islandRoots={};w.islandCanSleep={};w.islandMinTime={}
    w.joints={};w.jointNoCollide={};w.nextJoint=1
    return w
end

M._constants = table.freeze({ sphere=SHAPE_SPHERE, box=SHAPE_BOX, capsule=SHAPE_CAPSULE })
return table.freeze(M)
