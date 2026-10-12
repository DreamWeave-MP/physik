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
local SHAPE_MESH = 4
local SHAPE_HEIGHTFIELD = 5
local SHAPE_COMPOUND = 6
local JOINT_DISTANCE = 1
local JOINT_BALL = 2
local JOINT_HINGE = 3

local EVENT_ENTER = 1
local EVENT_STAY = 2
local EVENT_EXIT = 3

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
-- A byte for each residue of a packed cell: set while a cell with that residue is occupied, so
-- a clear byte says a cell is empty without a lookup in the table of cells.
local GridFilterSize = 65521
local GRID_KEY_DELTA = {}
for i=1,#GRID_DX do
    GRID_KEY_DELTA[i]=GRID_DX[i]*CELL_BASE2+GRID_DY[i]*CELL_BASE+GRID_DZ[i]
end

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


local function qMul(av,aw,bv,bw)
    return av*bw+bv*aw+vector.cross(av,bv),aw*bw-vector.dot(av,bv)
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


local function readVec3(buf, index)
    local o=index*12
    return vector.create(buffer.readf32(buf,o),buffer.readf32(buf,o+4),buffer.readf32(buf,o+8))
end

local function writeVec3(buf,index,v)
    local o=index*12
    buffer.writef32(buf,o,v.X); buffer.writef32(buf,o+4,v.Y); buffer.writef32(buf,o+8,v.Z)
end

local function buildMeshGrid(vbuf,ibuf,triCount,cellSize,maxCellsPerTriangle)
    local heads,nextNode,nodeTri,overflow={},{},{},{}
    local nodeCount=0
    for t=0,triCount-1 do
        local io=t*12
        local ia,ib,ic=buffer.readu32(ibuf,io),buffer.readu32(ibuf,io+4),buffer.readu32(ibuf,io+8)
        local a,b,c=readVec3(vbuf,ia),readVec3(vbuf,ib),readVec3(vbuf,ic)
        local mn=vector.create(math.min(a.X,b.X,c.X),math.min(a.Y,b.Y,c.Y),math.min(a.Z,b.Z,c.Z))
        local mx=vector.create(math.max(a.X,b.X,c.X),math.max(a.Y,b.Y,c.Y),math.max(a.Z,b.Z,c.Z))
        local x0,y0,z0=math.floor(mn.X/cellSize),math.floor(mn.Y/cellSize),math.floor(mn.Z/cellSize)
        local x1,y1,z1=math.floor(mx.X/cellSize),math.floor(mx.Y/cellSize),math.floor(mx.Z/cellSize)
        local cells=(x1-x0+1)*(y1-y0+1)*(z1-z0+1)
        local valid=cells<=maxCellsPerTriangle
        if valid then
            for x=x0,x1 do for y=y0,y1 do for z=z0,z1 do
                local key=packCell(x,y,z)
                if key==nil then valid=false; break end
                nodeCount+=1
                nodeTri[nodeCount]=t
                nextNode[nodeCount]=heads[key] or 0
                heads[key]=nodeCount
            end if not valid then break end end if not valid then break end end
        end
        if not valid then overflow[#overflow+1]=t end
    end
    return heads,nextNode,nodeTri,overflow
end

function M.triangleMesh(vertices,indices,opts)
    opts=opts or {}
    local vbuf,vcount
    if typeof(vertices)=="buffer" then
        assert(buffer.len(vertices)%12==0,"physik.triangleMesh: vertex buffer must contain packed f32 xyz triples")
        vbuf=vertices; vcount=buffer.len(vertices)/12
    else
        assert(type(vertices)=="table" and #vertices>=3,"physik.triangleMesh: vertices must be a vector array or packed f32 buffer")
        vcount=#vertices; vbuf=buffer.create(vcount*12)
        for i=1,vcount do assert(typeof(vertices[i])=="vector","physik.triangleMesh: every vertex must be a vector"); writeVec3(vbuf,i-1,vertices[i]) end
    end
    local ibuf,icount
    if indices==nil then
        assert(vcount%3==0,"physik.triangleMesh: without indices, vertex count must be divisible by 3")
        icount=vcount; ibuf=buffer.create(icount*4)
        for i=0,icount-1 do buffer.writeu32(ibuf,i*4,i) end
    elseif typeof(indices)=="buffer" then
        assert(buffer.len(indices)%12==0,"physik.triangleMesh: index buffer length must be divisible by 3 u32 indices")
        ibuf=indices; icount=buffer.len(indices)/4
    else
        assert(type(indices)=="table" and #indices%3==0,"physik.triangleMesh: indices must be a triangle-list table or u32 buffer")
        icount=#indices; ibuf=buffer.create(icount*4)
        for i=1,icount do
            local idx=indices[i]
            assert(type(idx)=="number" and idx%1==0 and idx>=1 and idx<=vcount,"physik.triangleMesh: table indices are 1-based and must reference a vertex")
            buffer.writeu32(ibuf,(i-1)*4,idx-1)
        end
    end
    local triCount=icount/3
    assert(triCount>=1,"physik.triangleMesh: mesh must contain at least one triangle")
    local mn=vector.create(math.huge,math.huge,math.huge)
    local mx=vector.create(-math.huge,-math.huge,-math.huge)
    for i=0,vcount-1 do
        local v=readVec3(vbuf,i)
        mn=vector.create(math.min(mn.X,v.X),math.min(mn.Y,v.Y),math.min(mn.Z,v.Z))
        mx=vector.create(math.max(mx.X,v.X),math.max(mx.Y,v.Y),math.max(mx.Z,v.Z))
    end
    local extent=mx-mn
    local longest=math.max(extent.X,extent.Y,extent.Z)
    local cellSize=opts.cellSize or math.max(longest/math.max(math.sqrt(triCount)*0.5,1),0.25)
    assert(type(cellSize)=="number" and cellSize>0,"physik.triangleMesh: cellSize must be > 0")
    local maxCells=opts.maxCellsPerTriangle or 64
    local heads,nextNode,nodeTri,overflow=buildMeshGrid(vbuf,ibuf,triCount,cellSize,maxCells)
    return table.freeze({kind=SHAPE_MESH,vertices=vbuf,indices=ibuf,vertexCount=vcount,triangleCount=triCount,localMin=mn,localMax=mx,cellSize=cellSize,heads=heads,nextNode=nextNode,nodeTri=nodeTri,overflow=overflow})
end

function M.heightfield(width,depth,heights,opts)
    opts=opts or {}
    assert(type(width)=="number" and width%1==0 and width>=2,"physik.heightfield: width must be an integer >= 2")
    assert(type(depth)=="number" and depth%1==0 and depth>=2,"physik.heightfield: depth must be an integer >= 2")
    local count=width*depth
    local hbuf
    if typeof(heights)=="buffer" then
        assert(buffer.len(heights)==count*4,"physik.heightfield: height buffer must contain width*depth f32 values")
        hbuf=heights
    else
        assert(type(heights)=="table" and #heights==count,"physik.heightfield: heights must contain width*depth values")
        hbuf=buffer.create(count*4)
        for i=1,count do assert(type(heights[i])=="number","physik.heightfield: heights must be numbers");buffer.writef32(hbuf,(i-1)*4,heights[i]) end
    end
    local sx=opts.cellSizeX or opts.cellSize or 1
    local sz=opts.cellSizeZ or opts.cellSize or 1
    assert(sx>0 and sz>0,"physik.heightfield: cell sizes must be > 0")
    local origin=opts.origin or ZERO
    assert(typeof(origin)=="vector","physik.heightfield: origin must be a vector")
    local minY,maxY=math.huge,-math.huge
    for i=0,count-1 do local h=buffer.readf32(hbuf,i*4);if h<minY then minY=h end;if h>maxY then maxY=h end end
    local mn=vector.create(origin.X,origin.Y+minY,origin.Z)
    local mx=vector.create(origin.X+(width-1)*sx,origin.Y+maxY,origin.Z+(depth-1)*sz)
    return table.freeze({kind=SHAPE_HEIGHTFIELD,width=width,depth=depth,heights=hbuf,cellSizeX=sx,cellSizeZ=sz,origin=origin,localMin=mn,localMax=mx})
end


local function primitiveLocalBounds(shape,p,qv,qw)
    if shape.kind==SHAPE_SPHERE then
        local e=vector.create(shape.radius,shape.radius,shape.radius)
        return p-e,p+e
    elseif shape.kind==SHAPE_BOX then
        local a0,a1,a2=qAxes(qv,qw)
        local h=shape.halfExtents
        local e=vector.create(math.abs(a0.X)*h.X+math.abs(a1.X)*h.Y+math.abs(a2.X)*h.Z,math.abs(a0.Y)*h.X+math.abs(a1.Y)*h.Y+math.abs(a2.Y)*h.Z,math.abs(a0.Z)*h.X+math.abs(a1.Z)*h.Y+math.abs(a2.Z)*h.Z)
        return p-e,p+e
    else
        local _,axis,_=qAxes(qv,qw)
        local e=vector.create(math.abs(axis.X)*shape.halfHeight+shape.radius,math.abs(axis.Y)*shape.halfHeight+shape.radius,math.abs(axis.Z)*shape.halfHeight+shape.radius)
        return p-e,p+e
    end
end

function M.compound(children)
    assert(type(children)=="table" and #children>=1,"physik.compound: children must be a non-empty array")
    local shapes,positions,qvs,qws={},{},{},{}
    local mn=vector.create(math.huge,math.huge,math.huge)
    local mx=vector.create(-math.huge,-math.huge,-math.huge)
    for i=1,#children do
        local child=children[i]
        assert(type(child)=="table" and type(child.shape)=="table","physik.compound: every child needs a shape")
        local shape=child.shape
        assert(shape.kind==SHAPE_SPHERE or shape.kind==SHAPE_BOX or shape.kind==SHAPE_CAPSULE,"physik.compound: children must be sphere, box, or capsule shapes")
        local p=child.position or ZERO
        assert(typeof(p)=="vector","physik.compound: child position must be a vector")
        local qv,qw=ZERO,1
        if child.rotation then
            assert(type(child.rotation)=="table" and typeof(child.rotation.v)=="vector" and type(child.rotation.w)=="number","physik.compound: child rotation must come from physik.quat/axisAngle")
            qv,qw=qNormalize(child.rotation.v,child.rotation.w)
        end
        shapes[i],positions[i],qvs[i],qws[i]=shape,p,qv,qw
        local cmn,cmx=primitiveLocalBounds(shape,p,qv,qw)
        mn=vector.create(math.min(mn.X,cmn.X),math.min(mn.Y,cmn.Y),math.min(mn.Z,cmn.Z))
        mx=vector.create(math.max(mx.X,cmx.X),math.max(mx.Y,cmx.Y),math.max(mx.Z,cmx.Z))
    end
    return table.freeze({kind=SHAPE_COMPOUND,count=#children,shapes=shapes,positions=positions,qv=qvs,qw=qws,localMin=mn,localMax=mx})
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
    elseif shape.kind == SHAPE_CAPSULE then
        -- Capsule: compact approximation by a box enclosing the swept sphere.
        local r, hh = shape.radius, shape.halfHeight
        local hx, hy, hz = r, hh + r, r
        local ix = (mass / 3) * (hy*hy + hz*hz)
        local iy = (mass / 3) * (hx*hx + hz*hz)
        local iz = (mass / 3) * (hx*hx + hy*hy)
        return vector.create(if ix > EPS then 1/ix else 0, if iy > EPS then 1/iy else 0, if iz > EPS then 1/iz else 0)
    elseif shape.kind == SHAPE_COMPOUND then
        -- Conservative compound inertia: enclosing cuboid about the declared body origin,
        -- including a parallel-axis term when the child bounds are offset.
        local h=(shape.localMax-shape.localMin)*0.5
        local c=(shape.localMax+shape.localMin)*0.5
        local ix=(mass/3)*(h.Y*h.Y+h.Z*h.Z)+mass*(c.Y*c.Y+c.Z*c.Z)
        local iy=(mass/3)*(h.X*h.X+h.Z*h.Z)+mass*(c.X*c.X+c.Z*c.Z)
        local iz=(mass/3)*(h.X*h.X+h.Y*h.Y)+mass*(c.X*c.X+c.Y*c.Y)
        return vector.create(if ix>EPS then 1/ix else 0,if iy>EPS then 1/iy else 0,if iz>EPS then 1/iz else 0)
    end
    return ZERO
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
    local p,dp,lo,hi,inv,t1,t2

    p,dp,lo,hi=a.X,d.X,-h.X,h.X
    if math.abs(dp)<=EPS then
        if p<lo or p>hi then return false end
    else
        inv=1/dp;t1=(lo-p)*inv;t2=(hi-p)*inv
        if t1>t2 then t1,t2=t2,t1 end
        if t1>tmin then tmin=t1 end
        if t2<tmax then tmax=t2 end
        if tmin>tmax then return false end
    end

    p,dp,lo,hi=a.Y,d.Y,-h.Y,h.Y
    if math.abs(dp)<=EPS then
        if p<lo or p>hi then return false end
    else
        inv=1/dp;t1=(lo-p)*inv;t2=(hi-p)*inv
        if t1>t2 then t1,t2=t2,t1 end
        if t1>tmin then tmin=t1 end
        if t2<tmax then tmax=t2 end
        if tmin>tmax then return false end
    end

    p,dp,lo,hi=a.Z,d.Z,-h.Z,h.Z
    if math.abs(dp)<=EPS then
        if p<lo or p>hi then return false end
    else
        inv=1/dp;t1=(lo-p)*inv;t2=(hi-p)*inv
        if t1>t2 then t1,t2=t2,t1 end
        if t1>tmin then tmin=t1 end
        if t2<tmax then tmax=t2 end
        if tmin>tmax then return false end
    end
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


local function closestPointTriangle(p,a,b,c)
    local ab,ac=p and (b-a),c-a
    local ap=p-a
    local d1,d2=vector.dot(ab,ap),vector.dot(ac,ap)
    if d1<=0 and d2<=0 then return a end
    local bp=p-b
    local d3,d4=vector.dot(ab,bp),vector.dot(ac,bp)
    if d3>=0 and d4<=d3 then return b end
    local vc=d1*d4-d3*d2
    if vc<=0 and d1>=0 and d3<=0 then
        local v=d1/(d1-d3)
        return a+ab*v
    end
    local cp=p-c
    local d5,d6=vector.dot(ab,cp),vector.dot(ac,cp)
    if d6>=0 and d5<=d6 then return c end
    local vb=d5*d2-d1*d6
    if vb<=0 and d2>=0 and d6<=0 then
        local w=d2/(d2-d6)
        return a+ac*w
    end
    local va=d3*d6-d5*d4
    if va<=0 and (d4-d3)>=0 and (d5-d6)>=0 then
        local w=(d4-d3)/((d4-d3)+(d5-d6))
        return b+(c-b)*w
    end
    local denom=1/(va+vb+vc)
    local v=vb*denom
    local w=vc*denom
    return a+ab*v+ac*w
end

local function pointInTriangle(p,a,b,c,n)
    local c0=vector.cross(b-a,p-a)
    local c1=vector.cross(c-b,p-b)
    local c2=vector.cross(a-c,p-c)
    local d0,d1,d2=vector.dot(c0,n),vector.dot(c1,n),vector.dot(c2,n)
    return (d0>=-1e-7 and d1>=-1e-7 and d2>=-1e-7) or (d0<=1e-7 and d1<=1e-7 and d2<=1e-7)
end

local function closestSegmentTriangle(s0,s1,a,b,c)
    local n=vector.cross(b-a,c-a)
    local n2=vector.dot(n,n)
    if n2>EPS then
        local d=s1-s0
        local denom=vector.dot(n,d)
        if math.abs(denom)>EPS then
            local t=vector.dot(n,a-s0)/denom
            if t>=0 and t<=1 then
                local p=s0+d*t
                if pointInTriangle(p,a,b,c,n) then return p,p end
            end
        end
    end
    local bestS=s0
    local bestT=closestPointTriangle(s0,a,b,c)
    local bestD=vector.dot(bestT-bestS,bestT-bestS)
    local q=closestPointTriangle(s1,a,b,c)
    local d2=vector.dot(q-s1,q-s1)
    if d2<bestD then bestD,bestS,bestT=d2,s1,q end
    local p0,p1=closestSegmentSegment(s0,s1,a,b)
    d2=vector.dot(p1-p0,p1-p0); if d2<bestD then bestD,bestS,bestT=d2,p0,p1 end
    p0,p1=closestSegmentSegment(s0,s1,b,c)
    d2=vector.dot(p1-p0,p1-p0); if d2<bestD then bestD,bestS,bestT=d2,p0,p1 end
    p0,p1=closestSegmentSegment(s0,s1,c,a)
    d2=vector.dot(p1-p0,p1-p0); if d2<bestD then bestS,bestT=p0,p1 end
    return bestS,bestT
end

local function sphereTriangle(sp,r,a,b,c,margin)
    local q=closestPointTriangle(sp,a,b,c)
    local d=q-sp
    local d2=vector.dot(d,d)
    local detect=r+(margin or 0)
    if d2>detect*detect then return false end
    local dist=math.sqrt(math.max(d2,EPS))
    local n
    if d2>EPS then n=d/dist else
        local tn=safeUnit(vector.cross(b-a,c-a),YAXIS)
        local center=(a+b+c)/3
        n=if vector.dot(center-sp,tn)>=0 then tn else -tn
    end
    return true,n,q,r-dist
end

local function capsuleTriangle(cp,cqv,cqw,r,hh,a,b,c,margin)
    local s0,s1=capsuleEndpoints(cp,cqv,cqw,hh)
    local ps,pt=closestSegmentTriangle(s0,s1,a,b,c)
    local d=pt-ps
    local d2=vector.dot(d,d)
    local detect=r+(margin or 0)
    if d2>detect*detect then return false end
    local dist=math.sqrt(math.max(d2,EPS))
    local n
    if d2>EPS then n=d/dist else
        local tn=safeUnit(vector.cross(b-a,c-a),YAXIS)
        local center=(a+b+c)/3
        n=if vector.dot(center-cp,tn)>=0 then tn else -tn
    end
    return true,n,pt,r-dist
end

local function triBoxAxis(axis,a,b,c,h,bestPen,bestAxis,margin)
    local l2=vector.dot(axis,axis)
    if l2<=1e-10 then return true,bestPen,bestAxis end
    local n=if math.abs(l2-1)<1e-5 then axis else axis/math.sqrt(l2)
    local pa,pb,pc=vector.dot(a,n),vector.dot(b,n),vector.dot(c,n)
    local tmin,tmax=math.min(pa,pb,pc),math.max(pa,pb,pc)
    local r=h.X*math.abs(n.X)+h.Y*math.abs(n.Y)+h.Z*math.abs(n.Z)
    local overlap=math.min(r,tmax)-math.max(-r,tmin)
    if overlap<=-(margin or 0) then return false,bestPen,bestAxis end
    if overlap<bestPen then
        bestPen=overlap
        local center=(tmin+tmax)*0.5
        bestAxis=if center>=0 then n else -n
    end
    return true,bestPen,bestAxis
end

local function boxTriangle(bp,bqv,bqw,h,wa,wb,wc,margin)
    local a=qInvRotate(bqv,bqw,wa-bp)
    local b=qInvRotate(bqv,bqw,wb-bp)
    local c=qInvRotate(bqv,bqw,wc-bp)
    local e0,e1,e2=b-a,c-b,a-c
    local bestPen,bestAxis=math.huge,nil
    local ok
    ok,bestPen,bestAxis=triBoxAxis(XAXIS,a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    ok,bestPen,bestAxis=triBoxAxis(YAXIS,a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    ok,bestPen,bestAxis=triBoxAxis(ZAXIS,a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    ok,bestPen,bestAxis=triBoxAxis(vector.cross(e0,e1),a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    ok,bestPen,bestAxis=triBoxAxis(vector.cross(e0,XAXIS),a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    ok,bestPen,bestAxis=triBoxAxis(vector.cross(e0,YAXIS),a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    ok,bestPen,bestAxis=triBoxAxis(vector.cross(e0,ZAXIS),a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    ok,bestPen,bestAxis=triBoxAxis(vector.cross(e1,XAXIS),a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    ok,bestPen,bestAxis=triBoxAxis(vector.cross(e1,YAXIS),a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    ok,bestPen,bestAxis=triBoxAxis(vector.cross(e1,ZAXIS),a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    ok,bestPen,bestAxis=triBoxAxis(vector.cross(e2,XAXIS),a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    ok,bestPen,bestAxis=triBoxAxis(vector.cross(e2,YAXIS),a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    ok,bestPen,bestAxis=triBoxAxis(vector.cross(e2,ZAXIS),a,b,c,h,bestPen,bestAxis,margin);if not ok then return false end
    local localPoint=closestPointTriangle(ZERO,a,b,c)
    -- A triangle has zero thickness, so the minimum SAT interval overlap can be nearly zero
    -- even while a box has crossed the triangle plane. SAT decides overlap; the triangle
    -- plane supplies the physically useful surface normal and support penetration depth.
    local triN=safeUnit(vector.cross(b-a,c-a),bestAxis or YAXIS)
    local planeD=(vector.dot(a,triN)+vector.dot(b,triN)+vector.dot(c,triN))/3
    if planeD<0 then triN=-triN;planeD=-planeD end
    local radius=h.X*math.abs(triN.X)+h.Y*math.abs(triN.Y)+h.Z*math.abs(triN.Z)
    local planePen=radius-planeD
    return true,qRotate(bqv,bqw,triN),bp+qRotate(bqv,bqw,localPoint),planePen
end

local function meshTriangle(shape,t)
    local io=t*12
    local ia,ib,ic=buffer.readu32(shape.indices,io),buffer.readu32(shape.indices,io+4),buffer.readu32(shape.indices,io+8)
    return readVec3(shape.vertices,ia),readVec3(shape.vertices,ib),readVec3(shape.vertices,ic)
end

local function heightVertex(shape,x,z)
    local h=buffer.readf32(shape.heights,(z*shape.width+x)*4)
    local o=shape.origin
    return vector.create(o.X+x*shape.cellSizeX,o.Y+h,o.Z+z*shape.cellSizeZ)
end

local function transformedBoundsLocal(qv,qw,p,mn,mx)
    local lmn=vector.create(math.huge,math.huge,math.huge)
    local lmx=vector.create(-math.huge,-math.huge,-math.huge)
    for xi=0,1 do for yi=0,1 do for zi=0,1 do
        local v=vector.create(if xi==0 then mn.X else mx.X,if yi==0 then mn.Y else mx.Y,if zi==0 then mn.Z else mx.Z)
        local l=qInvRotate(qv,qw,v-p)
        lmn=vector.create(math.min(lmn.X,l.X),math.min(lmn.Y,l.Y),math.min(lmn.Z,l.Z))
        lmx=vector.create(math.max(lmx.X,l.X),math.max(lmx.Y,l.Y),math.max(lmx.Z,l.Z))
    end end end
    return lmn,lmx
end


local function collidePrimitives(sa,ap,aqv,aqw,sb,bp,bqv,bqw,margin)
    local ok,n,p,pen
    if sa.kind==SHAPE_SPHERE and sb.kind==SHAPE_SPHERE then
        return sphereSphere(ap,sa.radius,bp,sb.radius,margin)
    elseif sa.kind==SHAPE_SPHERE and sb.kind==SHAPE_BOX then
        return sphereBox(ap,sa.radius,bp,bqv,bqw,sb.halfExtents,margin)
    elseif sa.kind==SHAPE_BOX and sb.kind==SHAPE_SPHERE then
        ok,n,p,pen=sphereBox(bp,sb.radius,ap,aqv,aqw,sa.halfExtents,margin);if ok then n=-n end;return ok,n,p,pen
    elseif sa.kind==SHAPE_BOX and sb.kind==SHAPE_BOX then
        local a0,a1,a2=qAxes(aqv,aqw);local b0,b1,b2=qAxes(bqv,bqw)
        ok,n,p,pen=boxBox(ap,a0,a1,a2,sa.halfExtents,bp,b0,b1,b2,sb.halfExtents,margin);return ok,n,p,pen
    elseif sa.kind==SHAPE_CAPSULE and sb.kind==SHAPE_SPHERE then
        return capsuleSphere(ap,aqv,aqw,sa.radius,sa.halfHeight,bp,sb.radius,margin)
    elseif sa.kind==SHAPE_SPHERE and sb.kind==SHAPE_CAPSULE then
        ok,n,p,pen=capsuleSphere(bp,bqv,bqw,sb.radius,sb.halfHeight,ap,sa.radius,margin);if ok then n=-n end;return ok,n,p,pen
    elseif sa.kind==SHAPE_CAPSULE and sb.kind==SHAPE_CAPSULE then
        return capsuleCapsule(ap,aqv,aqw,sa.radius,sa.halfHeight,bp,bqv,bqw,sb.radius,sb.halfHeight,margin)
    elseif sa.kind==SHAPE_CAPSULE and sb.kind==SHAPE_BOX then
        return capsuleBox(ap,aqv,aqw,sa.radius,sa.halfHeight,bp,bqv,bqw,sb.halfExtents,margin)
    elseif sa.kind==SHAPE_BOX and sb.kind==SHAPE_CAPSULE then
        ok,n,p,pen=capsuleBox(bp,bqv,bqw,sb.radius,sb.halfHeight,ap,aqv,aqw,sa.halfExtents,margin);if ok then n=-n end;return ok,n,p,pen
    end
    return false
end

local function primitiveTriangleContact(shape,p,qv,qw,a,b,c,margin)
    if shape.kind==SHAPE_SPHERE then return sphereTriangle(p,shape.radius,a,b,c,margin) end
    if shape.kind==SHAPE_CAPSULE then return capsuleTriangle(p,qv,qw,shape.radius,shape.halfHeight,a,b,c,margin) end
    if shape.kind==SHAPE_BOX then return boxTriangle(p,qv,qw,shape.halfExtents,a,b,c,margin) end
    return false
end

local function rayTriangleT(origin,dir,maxDist,a,b,c)
    local e1,e2=b-a,c-a
    local p=vector.cross(dir,e2)
    local det=vector.dot(e1,p)
    if math.abs(det)<1e-9 then return nil end
    local inv=1/det
    local tvec=origin-a
    local u=vector.dot(tvec,p)*inv
    if u<0 or u>1 then return nil end
    local q=vector.cross(tvec,e1)
    local v=vector.dot(dir,q)*inv
    if v<0 or u+v>1 then return nil end
    local t=vector.dot(e2,q)*inv
    if t<0 or t>maxDist then return nil end
    return t
end

local function triangleRayNormal(a,b,c,dir)
    local n=safeUnit(vector.cross(b-a,c-a),YAXIS)
    if vector.dot(n,dir)>0 then n=-n end
    return n
end

local function rayTriangle(origin,dir,maxDist,a,b,c)
    local t=rayTriangleT(origin,dir,maxDist,a,b,c)
    if not t then return nil end
    return t,triangleRayNormal(a,b,c,dir),origin+dir*t
end


local function sphereTriangleGap(startPos,motion,t,r2,a,b,c)
    local p=startPos+motion*t
    local q=closestPointTriangle(p,a,b,c)
    local d=p-q
    return vector.dot(d,d)-r2,p,q
end

local function sweepSphereTriangle(startPos,endPos,radius,a,b,c)
    local q0=closestPointTriangle(startPos,a,b,c)
    local d0=startPos-q0
    local r2=radius*radius
    if vector.dot(d0,d0)<=r2 then return nil end
    local motion=endPos-startPos
    local lo,hi=0,1
    for _=1,10 do
        local aT=(2*lo+hi)/3
        local bT=(lo+2*hi)/3
        local ga=sphereTriangleGap(startPos,motion,aT,r2,a,b,c)
        local gb=sphereTriangleGap(startPos,motion,bT,r2,a,b,c)
        if ga<gb then hi=bT else lo=aT end
    end
    local tMin=(lo+hi)*0.5
    if sphereTriangleGap(startPos,motion,tMin,r2,a,b,c)>0 then return nil end
    lo,hi=0,tMin
    for _=1,14 do
        local mid=(lo+hi)*0.5
        if sphereTriangleGap(startPos,motion,mid,r2,a,b,c)<=0 then hi=mid else lo=mid end
    end
    local t=hi
    local _,p,q=sphereTriangleGap(startPos,motion,t,r2,a,b,c)
    return t,safeUnit(p-q,safeUnit(vector.cross(b-a,c-a),YAXIS))
end

local function sphereSegmentGap(startPos,motion,t,r2,a,b)
    local p=startPos+motion*t
    local q=closestPointSegment(p,a,b)
    local d=p-q
    return vector.dot(d,d)-r2,p,q
end

local function sweepSphereSegment(startPos,endPos,radius,a,b)
    local q0=closestPointSegment(startPos,a,b)
    local d0=startPos-q0
    local r2=radius*radius
    if vector.dot(d0,d0)<=r2 then return nil end
    local motion=endPos-startPos
    local lo,hi=0,1
    for _=1,9 do
        local aT=(2*lo+hi)/3
        local bT=(lo+2*hi)/3
        if sphereSegmentGap(startPos,motion,aT,r2,a,b)<sphereSegmentGap(startPos,motion,bT,r2,a,b) then hi=bT else lo=aT end
    end
    local tMin=(lo+hi)*0.5
    if sphereSegmentGap(startPos,motion,tMin,r2,a,b)>0 then return nil end
    lo,hi=0,tMin
    for _=1,13 do local mid=(lo+hi)*0.5;if sphereSegmentGap(startPos,motion,mid,r2,a,b)<=0 then hi=mid else lo=mid end end
    local t=hi
    local _,p,q=sphereSegmentGap(startPos,motion,t,r2,a,b)
    return t,safeUnit(p-q,YAXIS)
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
    assert(type(shape) == "table" and (shape.kind == SHAPE_SPHERE or shape.kind == SHAPE_BOX or shape.kind == SHAPE_CAPSULE or shape.kind == SHAPE_MESH or shape.kind == SHAPE_HEIGHTFIELD or shape.kind == SHAPE_COMPOUND), "physik.addBody: invalid shape")
    local kind = bodyTypeCode(desc.type)
    assert(not ((shape.kind==SHAPE_MESH or shape.kind==SHAPE_HEIGHTFIELD) and kind==BODY_DYNAMIC), "physik.addBody: triangle meshes and heightfields must be static or kinematic")
    local mass = if kind == BODY_DYNAMIC then (desc.mass or 1) else 0
    assert(mass >= 0, "physik.addBody: mass must be >= 0")
    local i = self:_allocSlot()
    self.bodyType[i] = kind
    self.shape[i] = shape
    self.sphereRadius[i] = shape.kind == SHAPE_SPHERE and shape.radius
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
    self.sensor[i] = desc.sensor == true
    self.prevPos[i] = self.pos[i]
    self.sleeping[i] = false
    self.sleepTime[i] = 0
    self.aabbMin[i] = ZERO
    self.aabbMax[i] = ZERO
    self:_updateAabb(i)
    self.gridValid=false
    self.gridTopologyDirty=true
    if self.ccd[i] then self.ccdCount+=1 end
    return makeHandle(i,self.generation[i])
end

function World:removeBody(handle)
    local i = self:_check(handle)
    for id=1,self.nextJoint-1 do
        local j=self.joints[id]
        if j and (j.a==i or j.b==i) then self:_dropJoint(id,j) end
    end
    self.alive[i] = false
    self.aliveCount -= 1
    if self.ccd[i] then self.ccdCount-=1 end
    self.free[#self.free+1] = i
    self.userData[i] = nil
    self.sleeping[i] = false
    -- A destroyed slot must not turn into a "stay" event if the numeric slot is reused.
    for k=1,#self.contactPairPrev do
        local key=self.contactPairPrev[k]
        if self.contactPairA[key]==i or self.contactPairB[key]==i then
            self.contactPairActive[key]=nil
            self.contactPairA[key]=nil;self.contactPairB[key]=nil;self.contactPairSensor[key]=nil
        end
    end
    self.gridValid=false
    self.gridTopologyDirty=true
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
        sensor = self.sensor[i],
        userData = self.userData[i],
    }
end

function World:setPosition(handle, p, wake)
    local i = self:_check(handle)
    self.pos[i] = p
    self.gridValid=false
    if wake ~= false then self:wake(handle) end
end

function World:setRotation(handle,q,wake)
    local i=self:_check(handle)
    self.qv[i],self.qw[i]=qNormalize(q.v,q.w)
    self.axis0[i],self.axis1[i],self.axis2[i]=qAxes(self.qv[i],self.qw[i])
    self.gridValid=false
    if wake~=false then self:wake(handle) end
end


function World:driveKinematic(handle,targetPosition,targetRotation,dt)
    local i=self:_check(handle)
    assert(self.bodyType[i]==BODY_KINEMATIC,"physik.driveKinematic: body must be kinematic")
    assert(typeof(targetPosition)=="vector" and type(dt)=="number" and dt>0,"physik.driveKinematic: target position must be a vector and dt > 0")
    self.linVel[i]=(targetPosition-self.pos[i])/dt
    if targetRotation then
        assert(type(targetRotation)=="table" and typeof(targetRotation.v)=="vector" and type(targetRotation.w)=="number","physik.driveKinematic: rotation must come from physik.quat/axisAngle")
        local tv,tw=qNormalize(targetRotation.v,targetRotation.w)
        local dv,dw=qMul(tv,tw,-self.qv[i],self.qw[i])
        if dw<0 then dv,dw=-dv,-dw end
        dw=clamp(dw,-1,1)
        local angle=2*math.acos(dw)
        local s2=math.max(1-dw*dw,0)
        if angle>1e-7 and s2>1e-12 then self.angVel[i]=dv*(angle/(math.sqrt(s2)*dt)) else self.angVel[i]=ZERO end
    else
        self.angVel[i]=ZERO
    end
end

function World:shiftOrigin(delta)
    assert(typeof(delta)=="vector","physik.shiftOrigin: delta must be a vector")
    for i=1,self.high do
        if self.alive[i] then
            self.pos[i]-=delta
            self.prevPos[i]-=delta
            self.aabbMin[i]-=delta
            self.aabbMax[i]-=delta
        end
    end
    for c=1,self.contactCount do self.cp[c]-=delta end
    self.gridValid=false
    return self
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

function World:applyTorque(handle, torque)
    local i=self:_check(handle)
    if self.bodyType[i]~=BODY_DYNAMIC then return end
    self.torque[i]+=torque
    self.sleeping[i]=false
    self.sleepTime[i]=0
end

function World:applyAngularImpulse(handle, impulse)
    local i=self:_check(handle)
    if self.bodyType[i]~=BODY_DYNAMIC then return end
    self.angVel[i]+=applyInvInertia(self,i,impulse)
    self.sleeping[i]=false
    self.sleepTime[i]=0
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
    self.jointCount+=1
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
    self.jointCount+=1
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
    self.jointCount+=1
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
    self.jointCount-=1
end

function World:removeJoint(id)
    self:_dropJoint(id)
end

function World:_updateAabb(i)
    local shape = self.shape[i]
    local p = self.pos[i]
    if shape.kind == SHAPE_SPHERE then
        local r=shape.radius+self.contactMargin
        local e=vector.create(r,r,r)
        self.aabbMin[i],self.aabbMax[i]=p-e,p+e
        return
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
    elseif shape.kind == SHAPE_CAPSULE then
        local a0,axis,a2=qAxes(self.qv[i],self.qw[i])
        self.axis0[i],self.axis1[i],self.axis2[i]=a0,axis,a2
        local h = shape.halfHeight
        local r = shape.radius
        local e = vector.create(math.abs(axis.X)*h+r,math.abs(axis.Y)*h+r,math.abs(axis.Z)*h+r)
        self.aabbMin[i],self.aabbMax[i] = p-e,p+e
    else
        local a0,a1,a2=qAxes(self.qv[i],self.qw[i])
        self.axis0[i],self.axis1[i],self.axis2[i]=a0,a1,a2
        local lc=(shape.localMin+shape.localMax)*0.5
        local lh=(shape.localMax-shape.localMin)*0.5
        local wc=p+qRotate(self.qv[i],self.qw[i],lc)
        local e=vector.create(
            math.abs(a0.X)*lh.X+math.abs(a1.X)*lh.Y+math.abs(a2.X)*lh.Z,
            math.abs(a0.Y)*lh.X+math.abs(a1.Y)*lh.Y+math.abs(a2.Y)*lh.Z,
            math.abs(a0.Z)*lh.X+math.abs(a1.Z)*lh.Y+math.abs(a2.Z)*lh.Z
        )
        self.aabbMin[i],self.aabbMax[i]=wc-e,wc+e
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


function World:_pushGeomContact(n,p,pen)
    local count=self.geomContactCount
    for k=1,count do
        if vector.dot(self.geomPoint[k]-p,self.geomPoint[k]-p)<1e-6 then
            if pen>self.geomPen[k] then self.geomNormal[k],self.geomPoint[k],self.geomPen[k]=n,p,pen end
            return
        end
    end
    if count<4 then
        count+=1
        self.geomContactCount=count
        self.geomNormal[count],self.geomPoint[count],self.geomPen[count]=n,p,pen
        return
    end
    local shallow=1
    for k=2,4 do if self.geomPen[k]<self.geomPen[shallow] then shallow=k end end
    if pen>self.geomPen[shallow] then self.geomNormal[shallow],self.geomPoint[shallow],self.geomPen[shallow]=n,p,pen end
end

function World:_testGeomTriangle(prim,geom,la,lb,lc,flip)
    local gp,gqv,gqw=self.pos[geom],self.qv[geom],self.qw[geom]
    local a=gp+qRotate(gqv,gqw,la)
    local b=gp+qRotate(gqv,gqw,lb)
    local c=gp+qRotate(gqv,gqw,lc)
    local ps=self.shape[prim]
    if ps.kind==SHAPE_COMPOUND then
        local bp,bqv,bqw=self.pos[prim],self.qv[prim],self.qw[prim]
        for k=1,ps.count do
            local cv,cw=qMul(bqv,bqw,ps.qv[k],ps.qw[k])
            local cp=bp+qRotate(bqv,bqw,ps.positions[k])
            local ok,n,p,pen=primitiveTriangleContact(ps.shapes[k],cp,cv,cw,a,b,c,self.contactMargin)
            if ok then if flip then n=-n end;self:_pushGeomContact(n,p,pen) end
        end
        return
    end
    local ok,n,p,pen=primitiveTriangleContact(ps,self.pos[prim],self.qv[prim],self.qw[prim],a,b,c,self.contactMargin)
    if ok then if flip then n=-n end;self:_pushGeomContact(n,p,pen) end
end
function World:_meshTriangleOnce(prim,geom,shape,t,flip,epoch)
    if self.meshMarks[t]==epoch then return end
    self.meshMarks[t]=epoch
    local a,b,c=meshTriangle(shape,t)
    self:_testGeomTriangle(prim,geom,a,b,c,flip)
end

function World:_narrowMesh(prim,geom,flip)
    local shape=self.shape[geom]
    local lmn,lmx=transformedBoundsLocal(self.qv[geom],self.qw[geom],self.pos[geom],self.aabbMin[prim],self.aabbMax[prim])
    local cs=shape.cellSize
    local x0,y0,z0=math.floor(lmn.X/cs),math.floor(lmn.Y/cs),math.floor(lmn.Z/cs)
    local x1,y1,z1=math.floor(lmx.X/cs),math.floor(lmx.Y/cs),math.floor(lmx.Z/cs)
    local cells=(x1-x0+1)*(y1-y0+1)*(z1-z0+1)
    self.meshMarkEpoch+=1
    local epoch=self.meshMarkEpoch
    if cells>self.meshQueryCellLimit then
        for t=0,shape.triangleCount-1 do self:_meshTriangleOnce(prim,geom,shape,t,flip,epoch) end
        return
    end
    for x=x0,x1 do for y=y0,y1 do for z=z0,z1 do
        local key=packCell(x,y,z)
        if key then
            local node=shape.heads[key] or 0
            while node~=0 do
                self:_meshTriangleOnce(prim,geom,shape,shape.nodeTri[node],flip,epoch)
                node=shape.nextNode[node]
            end
        end
    end end end
    for k=1,#shape.overflow do self:_meshTriangleOnce(prim,geom,shape,shape.overflow[k],flip,epoch) end
end

function World:_narrowHeightfield(prim,geom,flip)
    local shape=self.shape[geom]
    local lmn,lmx=transformedBoundsLocal(self.qv[geom],self.qw[geom],self.pos[geom],self.aabbMin[prim],self.aabbMax[prim])
    local o=shape.origin
    local x0=clamp(math.floor((lmn.X-o.X)/shape.cellSizeX),0,shape.width-2)
    local x1=clamp(math.floor((lmx.X-o.X)/shape.cellSizeX),0,shape.width-2)
    local z0=clamp(math.floor((lmn.Z-o.Z)/shape.cellSizeZ),0,shape.depth-2)
    local z1=clamp(math.floor((lmx.Z-o.Z)/shape.cellSizeZ),0,shape.depth-2)
    if x1<x0 or z1<z0 then return end
    for z=z0,z1 do for x=x0,x1 do
        local v00=heightVertex(shape,x,z)
        local v10=heightVertex(shape,x+1,z)
        local v01=heightVertex(shape,x,z+1)
        local v11=heightVertex(shape,x+1,z+1)
        self:_testGeomTriangle(prim,geom,v00,v10,v11,flip)
        self:_testGeomTriangle(prim,geom,v00,v11,v01,flip)
    end end
end

function World:_narrowGeometry(a,b)
    local sa,sb=self.shape[a],self.shape[b]
    local aGeom=sa.kind==SHAPE_MESH or sa.kind==SHAPE_HEIGHTFIELD
    local bGeom=sb.kind==SHAPE_MESH or sb.kind==SHAPE_HEIGHTFIELD
    if aGeom and bGeom then return end
    local prim,geom,flip
    if aGeom then prim,geom,flip=b,a,true else prim,geom,flip=a,b,false end
    self.geomContactCount=0
    if self.shape[geom].kind==SHAPE_MESH then self:_narrowMesh(prim,geom,flip) else self:_narrowHeightfield(prim,geom,flip) end
    for k=1,self.geomContactCount do self:_addContact(a,b,self.geomNormal[k],self.geomPoint[k],self.geomPen[k]) end
end


function World:_testCompoundPrimitive(ashape,ap,aqv,aqw,bshape,bp,bqv,bqw)
    local ok,n,p,pen=collidePrimitives(ashape,ap,aqv,aqw,bshape,bp,bqv,bqw,self.contactMargin)
    if ok then self:_pushGeomContact(n,p,pen) end
end

function World:_narrowCompound(a,b)
    local sa,sb=self.shape[a],self.shape[b]
    local aComp=sa.kind==SHAPE_COMPOUND
    local bComp=sb.kind==SHAPE_COMPOUND
    self.geomContactCount=0
    if aComp and bComp then
        for i=1,sa.count do
            local aqv,aqw=qMul(self.qv[a],self.qw[a],sa.qv[i],sa.qw[i]);local ap=self.pos[a]+qRotate(self.qv[a],self.qw[a],sa.positions[i])
            for j=1,sb.count do
                local bqv,bqw=qMul(self.qv[b],self.qw[b],sb.qv[j],sb.qw[j]);local bp=self.pos[b]+qRotate(self.qv[b],self.qw[b],sb.positions[j])
                self:_testCompoundPrimitive(sa.shapes[i],ap,aqv,aqw,sb.shapes[j],bp,bqv,bqw)
            end
        end
    elseif aComp then
        for i=1,sa.count do
            local aqv,aqw=qMul(self.qv[a],self.qw[a],sa.qv[i],sa.qw[i]);local ap=self.pos[a]+qRotate(self.qv[a],self.qw[a],sa.positions[i])
            self:_testCompoundPrimitive(sa.shapes[i],ap,aqv,aqw,sb,self.pos[b],self.qv[b],self.qw[b])
        end
    else
        for j=1,sb.count do
            local bqv,bqw=qMul(self.qv[b],self.qw[b],sb.qv[j],sb.qw[j]);local bp=self.pos[b]+qRotate(self.qv[b],self.qw[b],sb.positions[j])
            self:_testCompoundPrimitive(sa,self.pos[a],self.qv[a],self.qw[a],sb.shapes[j],bp,bqv,bqw)
        end
    end
    for k=1,self.geomContactCount do self:_addContact(a,b,self.geomNormal[k],self.geomPoint[k],self.geomPen[k]) end
end
function World:_narrow(a,b)
    local sa,sb=self.shape[a],self.shape[b]
    if sa.kind==SHAPE_MESH or sa.kind==SHAPE_HEIGHTFIELD or sb.kind==SHAPE_MESH or sb.kind==SHAPE_HEIGHTFIELD then
        self:_narrowGeometry(a,b)
        return
    end
    if sa.kind==SHAPE_COMPOUND or sb.kind==SHAPE_COMPOUND then
        self:_narrowCompound(a,b)
        return
    end
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

--- Offers a pair from a broadphase to the narrow phase, unless the pair is filtered out or its
--- bounds are apart. Most pairs a broadphase offers are apart, and nothing else is read for those:
--- a sphere is bounded by its centre and radius, anything else by its box.
local function tryPair(self, first: number, second: number)
  local radii, pos = self.sphereRadius, self.pos
  local firstRadius, secondRadius = radii[first], radii[second]

  if firstRadius and secondRadius then
    local reach = firstRadius + secondRadius + self.contactMargin * 2
    local between = pos[second] - pos[first]
    if math.abs(between.X) > reach or math.abs(between.Y) > reach or math.abs(between.Z) > reach then return end
  elseif firstRadius or secondRadius then
    local sphere, other = first, second
    if secondRadius then sphere, other = second, first end

    local centre, reach = pos[sphere], (firstRadius or secondRadius) + self.contactMargin
    local low, high = self.aabbMin[other], self.aabbMax[other]
    if centre.X + reach < low.X or centre.X - reach > high.X then return end
    if centre.Y + reach < low.Y or centre.Y - reach > high.Y then return end
    if centre.Z + reach < low.Z or centre.Z - reach > high.Z then return end
  else
    local low, high = self.aabbMin[first], self.aabbMax[first]
    local otherLow, otherHigh = self.aabbMin[second], self.aabbMax[second]
    if otherLow.X > high.X or otherHigh.X < low.X then return end
    if otherLow.Y > high.Y or otherHigh.Y < low.Y then return end
    if otherLow.Z > high.Z or otherHigh.Z < low.Z then return end
  end

  if first == second or not self.alive[first] or not self.alive[second] then return end

  -- A pair is always narrowed lower slot first, so a contact matches the one cached for it.
  if second < first then
    first, second = second, first
  end

  -- Every entry in the table is a joint's, so with no joint there is nothing to look up.
  if self.jointCount > 0 and self.jointNoCollide[bodyPairKey(first, second)] then return end
  if self.bodyType[first] == BODY_STATIC and self.bodyType[second] == BODY_STATIC then return end

  local mask, group = self.mask, self.group
  if bit32.band(mask[first], group[second]) == 0 or bit32.band(mask[second], group[first]) == 0 then return end

  self.candidateCount += 1
  local before = self.contactCount
  self:_narrow(first, second)
  if self.contactCount > before then self.manifoldCount += 1 end
end



function World:_collideSap()
  self:_sortOrder()
  local {order, alive, aabbMin, aabbMax} = self
  local count = #order

  for position = 1, count do
    local first = order[position]
    if not alive[first] then continue end

    local reach = aabbMax[first].X

    for after = position + 1, count do
      local second = order[after]
      if not alive[second] then continue end
      if aabbMin[second].X > reach then break end

      tryPair(self, first, second)
    end
  end
end

--- Lists the cells holding more than one body, and every pair of occupied neighbouring cells, in
--- the order the cells were first occupied.
function World:_rebuildGridPairs()
  local {gridHeads: heads, gridKeys: keys, gridNext: nextIn, gridOccupied: occupied} = self
  local {gridDenseKeys: dense, gridPairA: pairA, gridPairB: pairB} = self
  local {gridKeyX: columns, gridKeyY: rows, gridKeyZ: layers} = self
  local denseCount, pairCount = 0, 0

  for index = 1, #keys do
    local cell = keys[index]
    local head = heads[cell] or 0

    if head ~= 0 and (nextIn[head] or 0) ~= 0 then
      denseCount += 1
      dense[denseCount] = cell
    end

    local column, row, layer = columns[index], rows[index], layers[index]
    local interior = column < CELL_OFFSET - 1
      and row > -CELL_OFFSET
      and row < CELL_OFFSET - 1
      and layer > -CELL_OFFSET
      and layer < CELL_OFFSET - 1

    for neighbour = 2, 14 do
      local key = if interior
        then cell + GRID_KEY_DELTA[neighbour]
        else packCell(column + GRID_DX[neighbour], row + GRID_DY[neighbour], layer + GRID_DZ[neighbour])

      -- Most neighbours are empty, and the filter says so without a lookup.
      if not key or buffer.readu8(occupied, key % GridFilterSize) == 0 or not heads[key] then continue end

      pairCount += 1
      pairA[pairCount], pairB[pairCount] = cell, key
    end
  end

  for index = #dense, denseCount + 1, -1 do
    dense[index] = nil
  end

  for index = #pairA, pairCount + 1, -1 do
    pairA[index], pairB[index] = nil, nil
  end

  self.gridDenseCount = denseCount
  self.gridPairCount = pairCount
end

function World:_buildGridIndex(refreshAabbs)
  local {gridHeads: heads, gridKeys: keys, gridNext: nextIn, gridLarge: large, gridOccupied: occupied} = self
  local {alive, sphereRadius: radii, pos, aabbMin, aabbMax, gridKey, gridX, gridY, gridZ} = self
  local cellSize, margin = self.broadphaseCellSize, self.contactMargin

  -- Most coherent worlds keep bodies in the same broad-phase cells for many frames. Validate
  -- membership first; when topology is unchanged, retain bucket chains and occupied-neighbor
  -- pairs instead of rebuilding Lua hash tables and probing thirteen empty neighbors per cell.
  if self.gridTopologyReady and not self.gridTopologyDirty then
    local same = true

    for body = 1, self.high do
      if not alive[body] then continue end

      local radius = radii[body]
      local width, height, depth

      if radius then
        width = (radius + margin) * 2
        height, depth = width, width
      else
        if refreshAabbs then self:_updateAabb(body) end

        local low, high = aabbMin[body], aabbMax[body]
        width, height, depth = high.X - low.X, high.Y - low.Y, high.Z - low.Z
      end

      if gridKey[body] then
        if width > cellSize or height > cellSize or depth > cellSize then
          same = false
          continue
        end

        local centre = pos[body]

        if
          math.floor(centre.X / cellSize) ~= gridX[body]
          or math.floor(centre.Y / cellSize) ~= gridY[body]
          or math.floor(centre.Z / cellSize) ~= gridZ[body]
        then
          same = false
        end
      elseif width <= cellSize and height <= cellSize and depth <= cellSize then
        local centre = pos[body]
        local column, row, layer = math.floor(centre.X / cellSize), math.floor(centre.Y / cellSize), math.floor(centre.Z / cellSize)
        if packCell(column, row, layer) then same = false end
      end
    end

    if same then
      self.gridValid = true
      return false
    end
  end

  if refreshAabbs and (not self.gridTopologyReady or self.gridTopologyDirty) then
    for body = 1, self.high do
      if alive[body] and not radii[body] then self:_updateAabb(body) end
    end
  end

  -- Every head belongs to a key, so emptying both tables forgets every cell.
  table.clear(heads)
  table.clear(keys)
  table.clear(large)
  buffer.fill(occupied, 0, 0)

  local {gridKeyX: columns, gridKeyY: rows, gridKeyZ: layers} = self
  local smallCount, keyCount = 0, 0

  for body = 1, self.high do
    if not alive[body] then
      gridKey[body] = nil
      continue
    end

    local radius = radii[body]
    local width, height, depth

    if radius then
      width = (radius + margin) * 2
      height, depth = width, width
    else
      local low, high = aabbMin[body], aabbMax[body]
      width, height, depth = high.X - low.X, high.Y - low.Y, high.Z - low.Z
    end

    local centre = pos[body]
    local column, row, layer = math.floor(centre.X / cellSize), math.floor(centre.Y / cellSize), math.floor(centre.Z / cellSize)
    local key = width <= cellSize and height <= cellSize and depth <= cellSize and packCell(column, row, layer)

    if not key then
      gridKey[body] = nil
      table.insert(large, body)
      continue
    end

    smallCount += 1
    gridX[body], gridY[body], gridZ[body] = column, row, layer
    gridKey[body] = key

    local head = heads[key]

    if not head then
      keyCount += 1
      keys[keyCount] = key
      columns[keyCount], rows[keyCount], layers[keyCount] = column, row, layer
      buffer.writeu8(occupied, key % GridFilterSize, 1)
    end

    nextIn[body] = head or 0
    heads[key] = body
  end

  self.gridSmallCount = smallCount
  self.gridTopologyReady = true
  self.gridTopologyDirty = false
  self:_rebuildGridPairs()
  self.gridValid = true

  return true
end

function World:_ensureQueryGrid()
    if self.gridValid then return end
    self:_buildGridIndex(true)
end

function World:_collideGrid()
  self:_buildGridIndex(true)
  local {gridHeads: heads, gridNext: nextIn, gridDenseKeys: dense, gridPairA: pairA, gridPairB: pairB} = self

  -- Same-cell pairs only exist in cells cached as non-singleton during the most recent topology
  -- rebuild. Their body membership is stable until a center crosses a cell boundary.
  for index = 1, self.gridDenseCount do
    local first = heads[dense[index]] or 0

    while first ~= 0 do
      local second = nextIn[first] or 0

      while second ~= 0 do
        tryPair(self, first, second)
        second = nextIn[second] or 0
      end

      first = nextIn[first] or 0
    end
  end

  -- Cross-cell work walks only occupied neighboring cell pairs. Empty-neighbor hash probes are
  -- paid once per topology change rather than every simulation tick.
  for index = 1, self.gridPairCount do
    local first = heads[pairA[index]] or 0
    local head = heads[pairB[index]] or 0

    while first ~= 0 do
      local second = head

      while second ~= 0 do
        tryPair(self, first, second)
        second = nextIn[second] or 0
      end

      first = nextIn[first] or 0
    end
  end

  -- Bodies wider than one grid cell stay out of the hash grid. Pair anything involving
  -- them with a coherent X sweep instead of the old O(L^2) fallback. Small-small pairs
  -- are skipped here because the grid already emitted them.
  local {gridLarge: large, alive, gridKey} = self
  local largeCount = #large
  if largeCount == 0 then return end

  if largeCount <= self.directLargeThreshold then
    for index = 1, largeCount do
      local first = large[index]

      for second = 1, self.high do
        if alive[second] and second ~= first and (gridKey[second] or second > first) then tryPair(self, first, second) end
      end
    end

    return
  end

  self:_sortOrder()
  local {order, aabbMin, aabbMax} = self
  local count = #order

  for position = 1, count do
    local first = order[position]
    if not alive[first] then continue end

    local reach = aabbMax[first].X
    local firstLarge = not gridKey[first]

    for after = position + 1, count do
      local second = order[after]
      if not alive[second] then continue end
      if aabbMin[second].X > reach then break end
      if firstLarge or not gridKey[second] then tryPair(self, first, second) end
    end
  end
end

function World:_collide()
    self.contactCount=0
    self.candidateCount=0
    self.manifoldCount=0
    self.cacheEpoch+=1
    if self.broadphase=="sap" then
        for i=1,self.high do if self.alive[i] then self:_updateAabb(i) end end
        self:_collideSap()
    else
        self:_collideGrid()
    end
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
        if not self.sensor[a] and not self.sensor[b] and self.bodyType[a]==BODY_DYNAMIC and self.bodyType[b]==BODY_DYNAMIC then
            unionBodies(parent,size,a,b)
        end
    end
    for id=1,self.nextJoint-1 do
        local j=self.joints[id]
        if j then
            local a,b=j.a,j.b
            if self.alive[a] and self.alive[b] and self.generation[a]==j.ga and self.generation[b]==j.gb
                and self.bodyType[a]==BODY_DYNAMIC and self.bodyType[b]==BODY_DYNAMIC then
                unionBodies(parent,size,a,b)
            end
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

--- Seeds a contact with the impulses of the cached contact nearest its anchors, if the pair
--- touched last step and one of its cached contacts faces the same way.
local function matchCachedContact(self, contact: number, firstAnchor: vector, secondAnchor: vector, tangent1: vector, tangent2: vector)
  if not self.warmStarting or self.cbias[contact] > 0 then return end

  local first, second = self.ca[contact], self.cb[contact]
  local entry = self.contactCache[pairKey(first, second)]
  local epoch = self.cacheEpoch
  if not entry or entry[CACHE_EPOCH] ~= epoch - 1 then return end
  if entry[CACHE_GA] ~= self.generation[first] or entry[CACHE_GB] ~= self.generation[second] then return end

  if entry[CACHE_MATCH_EPOCH] ~= epoch then
    entry[CACHE_MATCH_EPOCH] = epoch
    entry[CACHE_USED_MASK] = 0
  end

  local best, nearest = 0, self.manifoldMatchDistance * self.manifoldMatchDistance * 2
  local normal = self.cn[contact]
  local used = entry[CACHE_USED_MASK] or 0

  for slot = 1, entry[CACHE_COUNT] or 0 do
    local record = CACHE_HEADER + (slot - 1) * CACHE_STRIDE
    if bit32.band(used, bit32.lshift(1, slot - 1)) ~= 0 or not (vector.dot(entry[record + 3], normal) > 0.85) then continue end

    local firstOffset, secondOffset = firstAnchor - entry[record + 1], secondAnchor - entry[record + 2]
    local distance = vector.dot(firstOffset, firstOffset) + vector.dot(secondOffset, secondOffset)

    if distance < nearest then
      nearest = distance
      best = slot
    end
  end

  if best == 0 then return end

  local record = CACHE_HEADER + (best - 1) * CACHE_STRIDE
  local factor = self.warmStartFactor
  local tangent = entry[record + 5] * factor

  entry[CACHE_USED_MASK] = bit32.bor(used, bit32.lshift(1, best - 1))
  self.cnormalImpulse[contact] = entry[record + 4] * factor
  self.ctangentImpulse1[contact] = vector.dot(tangent, tangent1)
  self.ctangentImpulse2[contact] = vector.dot(tangent, tangent2)
  self.warmStartedCount += 1
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

function World:_buildContactEvents()
    self.contactEventEpoch+=1
    local epoch=self.contactEventEpoch
    local seen=self.contactPairSeen
    local active=self.contactPairActive
    local pairA,pairB,pairSensor=self.contactPairA,self.contactPairB,self.contactPairSensor
    local prev,nextKeys=self.contactPairPrev,self.contactPairNext
    for k=#nextKeys,1,-1 do nextKeys[k]=nil end
    local eventCount=0
    for c=1,self.contactCount do
        local a,b=self.ca[c],self.cb[c]
        local key=bodyPairKey(a,b)
        if seen[key]~=epoch then
            seen[key]=epoch
            nextKeys[#nextKeys+1]=key
            pairA[key],pairB[key]=a,b
            local sensor=self.sensor[a] or self.sensor[b]
            pairSensor[key]=sensor
            eventCount+=1
            self.eventA[eventCount],self.eventB[eventCount]=a,b
            self.eventKind[eventCount]=if active[key] then EVENT_STAY else EVENT_ENTER
            self.eventSensor[eventCount]=sensor
            active[key]=true
        end
    end
    for k=1,#prev do
        local key=prev[k]
        if seen[key]~=epoch and active[key] then
            local a,b=pairA[key],pairB[key]
            if a and b then
                eventCount+=1
                self.eventA[eventCount],self.eventB[eventCount]=a,b
                self.eventKind[eventCount]=EVENT_EXIT
                self.eventSensor[eventCount]=pairSensor[key]==true
            end
            active[key]=nil;pairA[key]=nil;pairB[key]=nil;pairSensor[key]=nil
        end
        prev[k]=nil
    end
    self.contactPairPrev,self.contactPairNext=nextKeys,prev
    self.contactEventCount=eventCount
end

function World:_wakeContacts()
    local speed2=self.wakeSpeed*self.wakeSpeed
    local wakePen=self.wakePenetration
    for c=1,self.contactCount do
        local a,b=self.ca[c],self.cb[c]
        if self.sensor[a] or self.sensor[b] then continue end
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
  self.warmStartedCount = 0
  local {sensor, invMass, sleeping, pos, qv, qw, friction, restitution, linVel, angVel} = self
  local {activeContacts: active, cactive: actives, ca: bodiesA, cb: bodiesB, cp: points, cn: normals, cpen: depths} = self
  local {cra: firstArms, crb: secondArms, cla: firstAnchors, clb: secondAnchors} = self
  local {cnormalMass: normalMasses, cnormalAngA: firstNormalTurns, cnormalAngB: secondNormalTurns} = self
  local {ctangent1: tangents1, ctangent2: tangents2, ctangentMass1: tangentMasses1, ctangentMass2: tangentMasses2} = self
  local {ctangentAng1A: firstTangentTurns1, ctangentAng1B: secondTangentTurns1} = self
  local {ctangentAng2A: firstTangentTurns2, ctangentAng2B: secondTangentTurns2} = self
  local {cfriction: frictions, crest: restitutions, cbias: biases} = self
  local {cnormalImpulse: normalImpulses, ctangentImpulse1: tangentImpulses1, ctangentImpulse2: tangentImpulses2} = self
  local {baumgarte, maxDepenetrationVelocity} = self
  local biasSlop, pace = math.max(self.penetrationSlop, self.baumgarteSlop), math.max(dt, 1e-5)
  local activeCount = 0

  for contact = 1, self.contactCount do
    local first, second = bodiesA[contact], bodiesB[contact]
    local firstMoves, secondMoves = invMass[first] > 0 and not sleeping[first], invMass[second] > 0 and not sleeping[second]
    local solved = not (sensor[first] or sensor[second]) and (firstMoves or secondMoves)
    actives[contact] = solved

    if not solved then
      normalMasses[contact] = 0
      tangentMasses1[contact] = 0
      tangentMasses2[contact] = 0
      continue
    end

    activeCount += 1
    active[activeCount] = contact

    local point, normal = points[contact], normals[contact]
    local firstArm, secondArm = point - pos[first], point - pos[second]
    local firstAnchor = qInvRotate(qv[first], qw[first], firstArm)
    local secondAnchor = qInvRotate(qv[second], qw[second], secondArm)

    firstArms[contact], secondArms[contact] = firstArm, secondArm
    firstAnchors[contact], secondAnchors[contact] = firstAnchor, secondAnchor

    local firstInverse, firstTurn = contactAxisResponse(self, first, firstArm, normal)
    local secondInverse, secondTurn = contactAxisResponse(self, second, secondArm, normal)
    local inverse = firstInverse + secondInverse

    normalMasses[contact] = if inverse > EPS then 1 / inverse else 0
    firstNormalTurns[contact], secondNormalTurns[contact] = firstTurn, secondTurn

    local tangent1, tangent2 = contactTangents(normal)
    tangents1[contact], tangents2[contact] = tangent1, tangent2

    firstInverse, firstTurn = contactAxisResponse(self, first, firstArm, tangent1)
    secondInverse, secondTurn = contactAxisResponse(self, second, secondArm, tangent1)
    inverse = firstInverse + secondInverse
    tangentMasses1[contact] = if inverse > EPS then 1 / inverse else 0
    firstTangentTurns1[contact], secondTangentTurns1[contact] = firstTurn, secondTurn

    firstInverse, firstTurn = contactAxisResponse(self, first, firstArm, tangent2)
    secondInverse, secondTurn = contactAxisResponse(self, second, secondArm, tangent2)
    inverse = firstInverse + secondInverse
    tangentMasses2[contact] = if inverse > EPS then 1 / inverse else 0
    firstTangentTurns2[contact], secondTangentTurns2[contact] = firstTurn, secondTurn

    frictions[contact] = math.sqrt(friction[first] * friction[second])
    restitutions[contact] = math.max(restitution[first], restitution[second])
    biases[contact] = math.min(math.max(depths[contact] - biasSlop, 0) * baumgarte / pace, maxDepenetrationVelocity)
    matchCachedContact(self, contact, firstAnchor, secondAnchor, tangent1, tangent2)
  end

  for index = #active, activeCount + 1, -1 do
    active[index] = nil
  end

  self.activeContactCount = activeCount

  -- Seed this frame's solve with last frame's accumulated impulses. This is done after
  -- every contact mass is prepared, so a warm start never changes the effective mass basis.
  for index = 1, activeCount do
    local contact = active[index]
    local normalImpulse = normalImpulses[contact]
    local tangentImpulse1, tangentImpulse2 = tangentImpulses1[contact], tangentImpulses2[contact]
    if normalImpulse == 0 and tangentImpulse1 == 0 and tangentImpulse2 == 0 then continue end

    local first, second = bodiesA[contact], bodiesB[contact]
    local linear = normals[contact] * normalImpulse + tangents1[contact] * tangentImpulse1 + tangents2[contact] * tangentImpulse2

    if invMass[first] > 0 and not sleeping[first] then
      linVel[first] -= linear * invMass[first]
      angVel[first] -= firstNormalTurns[contact] * normalImpulse
        + firstTangentTurns1[contact] * tangentImpulse1
        + firstTangentTurns2[contact] * tangentImpulse2
    end

    if invMass[second] > 0 and not sleeping[second] then
      linVel[second] += linear * invMass[second]
      angVel[second] += secondNormalTurns[contact] * normalImpulse
        + secondTangentTurns1[contact] * tangentImpulse1
        + secondTangentTurns2[contact] * tangentImpulse2
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
  local {alive, sleeping, invMass, linVel, angVel, restitutionThreshold, activeContacts: active} = self
  local {ca: bodiesA, cb: bodiesB, cn: normals, cra: firstArms, crb: secondArms} = self
  local {cnormalMass: normalMasses, cnormalAngA: firstNormalTurns, cnormalAngB: secondNormalTurns} = self
  local {ctangent1: tangents1, ctangent2: tangents2, ctangentMass1: tangentMasses1, ctangentMass2: tangentMasses2} = self
  local {ctangentAng1A: firstTangentTurns1, ctangentAng1B: secondTangentTurns1} = self
  local {ctangentAng2A: firstTangentTurns2, ctangentAng2B: secondTangentTurns2} = self
  local {cfriction: frictions, crest: restitutions, cbias: biases} = self
  local {cnormalImpulse: normalImpulses, ctangentImpulse1: tangentImpulses1, ctangentImpulse2: tangentImpulses2} = self

  for iteration = 1, self.velocityIterations do
    for index = 1, self.activeContactCount do
      local contact = active[index]
      local first, second = bodiesA[contact], bodiesB[contact]
      if not alive[first] or not alive[second] then continue end

      local normal = normals[contact]
      local firstArm, secondArm = firstArms[contact], secondArms[contact]
      local firstMoves, secondMoves = invMass[first] > 0 and not sleeping[first], invMass[second] > 0 and not sleeping[second]
      local firstVelocity = if sleeping[first] then ZERO else linVel[first] + vector.cross(angVel[first], firstArm)
      local secondVelocity = if sleeping[second] then ZERO else linVel[second] + vector.cross(angVel[second], secondArm)
      local normalMass = normalMasses[contact]

      if normalMass > 0 then
        local approach = vector.dot(secondVelocity - firstVelocity, normal)
        local bounce = if iteration == 1 and approach < -restitutionThreshold then -restitutions[contact] * approach else 0
        local held = normalImpulses[contact]
        local impulse = math.max(held + (math.max(bounce, biases[contact]) - approach) * normalMass, 0)
        local change = impulse - held

        normalImpulses[contact] = impulse

        if change ~= 0 then
          if firstMoves then
            linVel[first] -= normal * (invMass[first] * change)
            angVel[first] -= firstNormalTurns[contact] * change
          end

          if secondMoves then
            linVel[second] += normal * (invMass[second] * change)
            angVel[second] += secondNormalTurns[contact] * change
          end
        end
      end

      local tangentMass1, tangentMass2 = tangentMasses1[contact], tangentMasses2[contact]
      if not (tangentMass1 > 0 or tangentMass2 > 0) then continue end

      firstVelocity = if sleeping[first] then ZERO else linVel[first] + vector.cross(angVel[first], firstArm)
      secondVelocity = if sleeping[second] then ZERO else linVel[second] + vector.cross(angVel[second], secondArm)

      local slip = secondVelocity - firstVelocity
      local tangent1, tangent2 = tangents1[contact], tangents2[contact]
      local held1, held2 = tangentImpulses1[contact], tangentImpulses2[contact]
      local impulse1, impulse2 = held1 - vector.dot(slip, tangent1) * tangentMass1, held2 - vector.dot(slip, tangent2) * tangentMass2
      local limit = frictions[contact] * normalImpulses[contact]
      local size = impulse1 * impulse1 + impulse2 * impulse2

      if size > limit * limit and size > EPS then
        local scale = limit / math.sqrt(size)
        impulse1, impulse2 = impulse1 * scale, impulse2 * scale
      end

      local change1, change2 = impulse1 - held1, impulse2 - held2
      tangentImpulses1[contact], tangentImpulses2[contact] = impulse1, impulse2
      if change1 == 0 and change2 == 0 then continue end

      local linear = tangent1 * change1 + tangent2 * change2

      if firstMoves then
        linVel[first] -= linear * invMass[first]
        angVel[first] -= firstTangentTurns1[contact] * change1 + firstTangentTurns2[contact] * change2
      end

      if secondMoves then
        linVel[second] += linear * invMass[second]
        angVel[second] += secondTangentTurns1[contact] * change1 + secondTangentTurns2[contact] * change2
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
    for id=1,self.nextJoint-1 do
        local j=self.joints[id]
        if not j then continue end
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


function World:_sweepSphereMesh(geom,startPos,endPos,radius)
    local shape=self.shape[geom]
    local gqv,gqw,gp=self.qv[geom],self.qw[geom],self.pos[geom]
    local ls=qInvRotate(gqv,gqw,startPos-gp)
    local le=qInvRotate(gqv,gqw,endPos-gp)
    local e=vector.create(radius,radius,radius)
    local mn=vector.create(math.min(ls.X,le.X),math.min(ls.Y,le.Y),math.min(ls.Z,le.Z))-e
    local mx=vector.create(math.max(ls.X,le.X),math.max(ls.Y,le.Y),math.max(ls.Z,le.Z))+e
    local cs=shape.cellSize
    local x0,y0,z0=math.floor(mn.X/cs),math.floor(mn.Y/cs),math.floor(mn.Z/cs)
    local x1,y1,z1=math.floor(mx.X/cs),math.floor(mx.Y/cs),math.floor(mx.Z/cs)
    local cells=(x1-x0+1)*(y1-y0+1)*(z1-z0+1)
    self.meshMarkEpoch+=1
    local epoch=self.meshMarkEpoch
    local bestT,bestN=1,nil
    if cells>self.meshQueryCellLimit then
        for t=0,shape.triangleCount-1 do
            self.meshMarks[t]=epoch
            local a,b,c=meshTriangle(shape,t)
            local tt,n=sweepSphereTriangle(ls,le,radius,a,b,c)
            if tt and tt<bestT then bestT,bestN=tt,n end
        end
    else
        for x=x0,x1 do for y=y0,y1 do for z=z0,z1 do
            local key=packCell(x,y,z)
            if key then
                local node=shape.heads[key] or 0
                while node~=0 do
                    local t=shape.nodeTri[node]
                    if self.meshMarks[t]~=epoch then
                        self.meshMarks[t]=epoch
                        local a,b,c=meshTriangle(shape,t)
                        local tt,n=sweepSphereTriangle(ls,le,radius,a,b,c)
                        if tt and tt<bestT then bestT,bestN=tt,n end
                    end
                    node=shape.nextNode[node]
                end
            end
        end end end
        for k=1,#shape.overflow do
            local t=shape.overflow[k]
            if self.meshMarks[t]~=epoch then
                self.meshMarks[t]=epoch
                local a,b,c=meshTriangle(shape,t)
                local tt,n=sweepSphereTriangle(ls,le,radius,a,b,c)
                if tt and tt<bestT then bestT,bestN=tt,n end
            end
        end
    end
    if not bestN then return nil end
    return bestT,qRotate(gqv,gqw,bestN)
end

function World:_sweepSphereHeightfield(geom,startPos,endPos,radius)
    local shape=self.shape[geom]
    local gqv,gqw,gp=self.qv[geom],self.qw[geom],self.pos[geom]
    local ls=qInvRotate(gqv,gqw,startPos-gp)
    local le=qInvRotate(gqv,gqw,endPos-gp)
    local o=shape.origin
    local x0=clamp(math.floor((math.min(ls.X,le.X)-radius-o.X)/shape.cellSizeX),0,shape.width-2)
    local x1=clamp(math.floor((math.max(ls.X,le.X)+radius-o.X)/shape.cellSizeX),0,shape.width-2)
    local z0=clamp(math.floor((math.min(ls.Z,le.Z)-radius-o.Z)/shape.cellSizeZ),0,shape.depth-2)
    local z1=clamp(math.floor((math.max(ls.Z,le.Z)+radius-o.Z)/shape.cellSizeZ),0,shape.depth-2)
    local bestT,bestN=1,nil
    for z=z0,z1 do for x=x0,x1 do
        local v00=heightVertex(shape,x,z)
        local v10=heightVertex(shape,x+1,z)
        local v01=heightVertex(shape,x,z+1)
        local v11=heightVertex(shape,x+1,z+1)
        local t,n=sweepSphereTriangle(ls,le,radius,v00,v10,v11)
        if t and t<bestT then bestT,bestN=t,n end
        t,n=sweepSphereTriangle(ls,le,radius,v00,v11,v01)
        if t and t<bestT then bestT,bestN=t,n end
    end end
    if not bestN then return nil end
    return bestT,qRotate(gqv,gqw,bestN)
end


function World:_sweepSphereCompound(geom,startPos,endPos,radius)
    local shape=self.shape[geom]
    local bp,bqv,bqw=self.pos[geom],self.qv[geom],self.qw[geom]
    local bestT,bestN=1,nil
    for k=1,shape.count do
        local cp=bp+qRotate(bqv,bqw,shape.positions[k])
        local cqv,cqw=qMul(bqv,bqw,shape.qv[k],shape.qw[k])
        local child=shape.shapes[k]
        local t,n
        if child.kind==SHAPE_SPHERE then
            t,n=sweepSphereSphere(startPos,endPos,radius,cp,child.radius)
        elseif child.kind==SHAPE_BOX then
            t,n=sweepSphereBox(startPos,endPos,radius,cp,cqv,cqw,child.halfExtents)
        else
            local a,b=capsuleEndpoints(cp,cqv,cqw,child.halfHeight)
            t,n=sweepSphereSegment(startPos,endPos,radius+child.radius,a,b)
        end
        if t and t<bestT then bestT,bestN=t,n end
    end
    if not bestN then return nil end
    return bestT,bestN
end

function World:_ccd()
    self.ccdHitCount=0
    if self.ccdCount==0 then return end
    for i=1,self.high do
        if not self.alive[i] or self.bodyType[i]~=BODY_DYNAMIC or self.sleeping[i] or not self.ccd[i] then continue end
        local shape=self.shape[i]
        if shape.kind~=SHAPE_SPHERE then continue end -- CCD is intentionally sphere-first.
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
            elseif target.kind==SHAPE_MESH then
                t,n=self:_sweepSphereMesh(j,start,finish,shape.radius)
            elseif target.kind==SHAPE_HEIGHTFIELD then
                t,n=self:_sweepSphereHeightfield(j,start,finish,shape.radius)
            elseif target.kind==SHAPE_COMPOUND then
                t,n=self:_sweepSphereCompound(j,start,finish,shape.radius)
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


function World:_sleepIndependent(dt)
    if not self.allowSleep then return end
    local lin2=self.sleepLinear*self.sleepLinear
    local ang2=self.sleepAngular*self.sleepAngular
    local count=0
    for i=1,self.high do
        if not self.alive[i] or self.bodyType[i]~=BODY_DYNAMIC then continue end
        count+=1
        if self.sleeping[i] then continue end
        if vector.dot(self.linVel[i],self.linVel[i])>=lin2 or vector.dot(self.angVel[i],self.angVel[i])>=ang2 then
            self.sleepTime[i]=0
        else
            local t=self.sleepTime[i]+dt
            self.sleepTime[i]=t
            if t>=self.sleepDelay then
                self.sleeping[i]=true
                self.linVel[i],self.angVel[i]=ZERO,ZERO
            end
        end
    end
    self.islandCount=count
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
        if self.trackContactEvents then self:_buildContactEvents() else self.contactEventCount=0 end
        if self.contactCount==0 and self.jointCount==0 then
            -- With no contact/joint graph every dynamic body is an independent island. Avoid
            -- building union-find state merely to rediscover that fact.
            self.activeContactCount=0
            self.warmStartedCount=0
            self:_sleepIndependent(h)
            if self.cacheEpoch%64==0 then
                for key,e in self.contactCache do
                    if self.cacheEpoch-e[CACHE_EPOCH]>2 then self.contactCache[key]=nil end
                end
            end
        else
            self:_wakeContacts()
            self:_buildIslands()
            self:_prepareContacts(h)
            self:_solveContacts(h)
            self:_storeContactCache()
            self:_correctPositions()
            self:_sleep(h)
        end
        -- Position correction can move bodies after the collision broadphase was built.
        -- Queries rebuild lazily only if a query is actually issued.
        self.gridValid=false
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

local function rayCenteredSlab(o,d,h,tmin,tmax,axis,hitAxis,hitSign)
    if math.abs(d)<EPS then return math.abs(o)<=h,tmin,tmax,hitAxis,hitSign end
    local inv=1/d
    local t1=(-h-o)*inv
    local t2=(h-o)*inv
    local sign=-1
    if t1>t2 then t1,t2=t2,t1;sign=1 end
    if t1>tmin then tmin=t1;hitAxis=axis;hitSign=sign end
    if t2<tmax then tmax=t2 end
    return tmin<=tmax,tmin,tmax,hitAxis,hitSign
end

local function rayBox(origin,dir,maxDist,center,qv,qw,half)
    local o=qInvRotate(qv,qw,origin-center)
    local d=qInvRotate(qv,qw,dir)
    local tmin,tmax=0,maxDist
    local hitAxis,hitSign=1,1
    local ok
    ok,tmin,tmax,hitAxis,hitSign=rayCenteredSlab(o.X,d.X,half.X,tmin,tmax,1,hitAxis,hitSign);if not ok then return nil end
    ok,tmin,tmax,hitAxis,hitSign=rayCenteredSlab(o.Y,d.Y,half.Y,tmin,tmax,2,hitAxis,hitSign);if not ok then return nil end
    ok,tmin,tmax,hitAxis,hitSign=rayCenteredSlab(o.Z,d.Z,half.Z,tmin,tmax,3,hitAxis,hitSign);if not ok then return nil end
    if tmin<0 or tmin>maxDist then return nil end
    local nl=if hitAxis==1 then vector.create(hitSign,0,0) elseif hitAxis==2 then vector.create(0,hitSign,0) else vector.create(0,0,hitSign)
    return tmin,qRotate(qv,qw,nl),origin+dir*tmin
end

local function rayRangeSlab(o,d,lo,hi,tmin,tmax)
    if math.abs(d)<EPS then return o>=lo and o<=hi,tmin,tmax end
    local inv=1/d
    local a,b=(lo-o)*inv,(hi-o)*inv
    if a>b then a,b=b,a end
    if a>tmin then tmin=a end
    if b<tmax then tmax=b end
    return tmin<=tmax,tmin,tmax
end

local function rayAabbInterval(origin,dir,maxDist,mn,mx)
    local tmin,tmax=0,maxDist
    local ok
    ok,tmin,tmax=rayRangeSlab(origin.X,dir.X,mn.X,mx.X,tmin,tmax);if not ok then return nil end
    ok,tmin,tmax=rayRangeSlab(origin.Y,dir.Y,mn.Y,mx.Y,tmin,tmax);if not ok then return nil end
    ok,tmin,tmax=rayRangeSlab(origin.Z,dir.Z,mn.Z,mx.Z,tmin,tmax);if not ok then return nil end
    return tmin,tmax
end

function World:_rayMeshTriangleT(shape,t,origin,dir,bestT,epoch)
    if self.meshMarks[t]==epoch then return nil end
    self.meshMarks[t]=epoch
    local a,b,c=meshTriangle(shape,t)
    return rayTriangleT(origin,dir,bestT,a,b,c)
end

function World:_rayMesh(i,origin,dir,maxDist)
    local shape=self.shape[i]
    local lo=qInvRotate(self.qv[i],self.qw[i],origin-self.pos[i])
    local ld=qInvRotate(self.qv[i],self.qw[i],dir)
    local enter,exit=rayAabbInterval(lo,ld,maxDist,shape.localMin,shape.localMax)
    if not enter then return nil end
    local a=lo+ld*math.max(enter,0)
    local b=lo+ld*exit
    local mn=vector.create(math.min(a.X,b.X),math.min(a.Y,b.Y),math.min(a.Z,b.Z))
    local mx=vector.create(math.max(a.X,b.X),math.max(a.Y,b.Y),math.max(a.Z,b.Z))
    local cs=shape.cellSize
    local x0,y0,z0=math.floor(mn.X/cs),math.floor(mn.Y/cs),math.floor(mn.Z/cs)
    local x1,y1,z1=math.floor(mx.X/cs),math.floor(mx.Y/cs),math.floor(mx.Z/cs)
    local cells=(x1-x0+1)*(y1-y0+1)*(z1-z0+1)
    self.meshMarkEpoch+=1
    local epoch=self.meshMarkEpoch
    local bestT,bestTri=maxDist,nil
    if cells>self.meshQueryCellLimit then
        for t=0,shape.triangleCount-1 do
            local tt=self:_rayMeshTriangleT(shape,t,lo,ld,bestT,epoch)
            if tt and tt<bestT then bestT,bestTri=tt,t end
        end
    else
        for x=x0,x1 do for y=y0,y1 do for z=z0,z1 do
            local key=packCell(x,y,z)
            if key then
                local node=shape.heads[key] or 0
                while node~=0 do
                    local tri=shape.nodeTri[node]
                    local tt=self:_rayMeshTriangleT(shape,tri,lo,ld,bestT,epoch)
                    if tt and tt<bestT then bestT,bestTri=tt,tri end
                    node=shape.nextNode[node]
                end
            end
        end end end
        for k=1,#shape.overflow do
            local tri=shape.overflow[k]
            local tt=self:_rayMeshTriangleT(shape,tri,lo,ld,bestT,epoch)
            if tt and tt<bestT then bestT,bestTri=tt,tri end
        end
    end
    if bestTri==nil then return nil end
    local ta,tb,tc=meshTriangle(shape,bestTri)
    local bestN=triangleRayNormal(ta,tb,tc,ld)
    return bestT,qRotate(self.qv[i],self.qw[i],bestN),origin+dir*bestT
end

function World:_rayHeightfield(i,origin,dir,maxDist)
    local shape=self.shape[i]
    local lo=qInvRotate(self.qv[i],self.qw[i],origin-self.pos[i])
    local ld=qInvRotate(self.qv[i],self.qw[i],dir)
    local enter,exit=rayAabbInterval(lo,ld,maxDist,shape.localMin,shape.localMax)
    if not enter then return nil end
    local a=lo+ld*math.max(enter,0)
    local b=lo+ld*exit
    local o=shape.origin
    local x0=clamp(math.floor((math.min(a.X,b.X)-o.X)/shape.cellSizeX),0,shape.width-2)
    local x1=clamp(math.floor((math.max(a.X,b.X)-o.X)/shape.cellSizeX),0,shape.width-2)
    local z0=clamp(math.floor((math.min(a.Z,b.Z)-o.Z)/shape.cellSizeZ),0,shape.depth-2)
    local z1=clamp(math.floor((math.max(a.Z,b.Z)-o.Z)/shape.cellSizeZ),0,shape.depth-2)
    local bestT,bestA,bestB,bestC=maxDist,nil,nil,nil
    for z=z0,z1 do for x=x0,x1 do
        local v00=heightVertex(shape,x,z)
        local v10=heightVertex(shape,x+1,z)
        local v01=heightVertex(shape,x,z+1)
        local v11=heightVertex(shape,x+1,z+1)
        local t=rayTriangleT(lo,ld,bestT,v00,v10,v11)
        if t and t<bestT then bestT,bestA,bestB,bestC=t,v00,v10,v11 end
        t=rayTriangleT(lo,ld,bestT,v00,v11,v01)
        if t and t<bestT then bestT,bestA,bestB,bestC=t,v00,v11,v01 end
    end end
    if not bestA then return nil end
    local bestN=triangleRayNormal(bestA,bestB,bestC,ld)
    return bestT,qRotate(self.qv[i],self.qw[i],bestN),origin+dir*bestT
end


local function rayPrimitiveShape(shape,p,qv,qw,origin,dir,maxDist)
    if shape.kind==SHAPE_SPHERE then return raySphere(origin,dir,maxDist,p,shape.radius) end
    if shape.kind==SHAPE_BOX then return rayBox(origin,dir,maxDist,p,qv,qw,shape.halfExtents) end
    local a,b=capsuleEndpoints(p,qv,qw,shape.halfHeight)
    local m=(a+b)*0.5
    local ta,na,pa=raySphere(origin,dir,maxDist,a,shape.radius)
    local tm,nm,pm=raySphere(origin,dir,maxDist,m,shape.radius)
    local tb,nb,pb=raySphere(origin,dir,maxDist,b,shape.radius)
    local t,n,pt=ta,na,pa
    if tm and (not t or tm<t) then t,n,pt=tm,nm,pm end
    if tb and (not t or tb<t) then t,n,pt=tb,nb,pb end
    return t,n,pt
end

function World:_rayCompound(i,origin,dir,maxDist)
    local shape=self.shape[i]
    local bp,bqv,bqw=self.pos[i],self.qv[i],self.qw[i]
    local bestT,bestN,bestP=maxDist,nil,nil
    for k=1,shape.count do
        local cp=bp+qRotate(bqv,bqw,shape.positions[k])
        local cqv,cqw=qMul(bqv,bqw,shape.qv[k],shape.qw[k])
        local t,n,p=rayPrimitiveShape(shape.shapes[k],cp,cqv,cqw,origin,dir,bestT)
        if t and t<bestT then bestT,bestN,bestP=t,n,p end
    end
    if not bestN then return nil end
    return bestT,bestN,bestP
end

function World:_raycastBody(i,origin,dir,bestT,mask)
    if not self.alive[i] or bit32.band(mask,self.group[i])==0 then return nil end
    local mn,mx=self.aabbMin[i],self.aabbMax[i]
    if not rayAabbInterval(origin,dir,bestT,mn,mx) then return nil end
    local s=self.shape[i]
    if s.kind==SHAPE_SPHERE then
        return raySphere(origin,dir,bestT,self.pos[i],s.radius)
    elseif s.kind==SHAPE_BOX then
        return rayBox(origin,dir,bestT,self.pos[i],self.qv[i],self.qw[i],s.halfExtents)
    elseif s.kind==SHAPE_MESH then
        return self:_rayMesh(i,origin,dir,bestT)
    elseif s.kind==SHAPE_HEIGHTFIELD then
        return self:_rayHeightfield(i,origin,dir,bestT)
    elseif s.kind==SHAPE_COMPOUND then
        return self:_rayCompound(i,origin,dir,bestT)
    end
    local a,b=capsuleEndpoints(self.pos[i],self.qv[i],self.qw[i],s.halfHeight)
    local m=(a+b)*0.5
    local ta,na,pa=raySphere(origin,dir,bestT,a,s.radius)
    local tm,nm,pm=raySphere(origin,dir,bestT,m,s.radius)
    local tb,nb,pb=raySphere(origin,dir,bestT,b,s.radius)
    local t,n,p=ta,na,pa
    if tm and (not t or tm<t) then t,n,p=tm,nm,pm end
    if tb and (not t or tb<t) then t,n,p=tb,nb,pb end
    return t,n,p
end

function World:_raycastRawScan(origin,dir,maxDistance,mask)
    local bestT,bestI,bestN,bestP=maxDistance,nil,nil,nil
    for i=1,self.high do
        local t,n,p=self:_raycastBody(i,origin,dir,bestT,mask)
        if t and t<bestT then bestT,bestI,bestN,bestP=t,i,n,p end
    end
    return bestI,bestT,bestP,bestN
end

function World:_raycastRaw(origin,direction,maxDistance,mask)
    maxDistance=maxDistance or math.huge
    mask=if mask==nil then 0xffffffff else mask
    local dl2=vector.dot(direction,direction)
    if dl2<=EPS then return nil end
    local dir=direction/math.sqrt(dl2)
    self:_ensureQueryGrid()

    -- Oversized bodies (world meshes, terrain, large compounds) are normally few. Test them first;
    -- a world hit gives the grid walk a finite early-out distance even for an unbounded ray.
    local bestT,bestI,bestN,bestP=maxDistance,nil,nil,nil
    local large=self.gridLarge
    for k=1,#large do
        local i=large[k]
        local t,n,p=self:_raycastBody(i,origin,dir,bestT,mask)
        if t and t<bestT then bestT,bestI,bestN,bestP=t,i,n,p end
    end
    if self.gridSmallCount==0 then return bestI,bestT,bestP,bestN end

    -- If the ray is effectively unbounded and oversized geometry did not cap it, a DDA could walk
    -- forever. Fall back to the linear scan in that unusual case.
    if bestT==math.huge then
        return self:_raycastRawScan(origin,dir,maxDistance,mask)
    end

    local cs=self.broadphaseCellSize
    local estimated=(math.abs(dir.X)+math.abs(dir.Y)+math.abs(dir.Z))*bestT/cs+3
    if estimated>self.queryRayCellLimit then
        return self:_raycastRawScan(origin,dir,maxDistance,mask)
    end

    self.queryMarkEpoch+=1
    local epoch=self.queryMarkEpoch
    local marks=self.queryMarks
    local heads,nextIn=self.gridHeads,self.gridNext
    local cx,cy,cz=math.floor(origin.X/cs),math.floor(origin.Y/cs),math.floor(origin.Z/cs)
    local sx=if dir.X>EPS then 1 elseif dir.X<-EPS then -1 else 0
    local sy=if dir.Y>EPS then 1 elseif dir.Y<-EPS then -1 else 0
    local sz=if dir.Z>EPS then 1 elseif dir.Z<-EPS then -1 else 0
    local inf=math.huge
    local tdx=if sx~=0 then cs/math.abs(dir.X) else inf
    local tdy=if sy~=0 then cs/math.abs(dir.Y) else inf
    local tdz=if sz~=0 then cs/math.abs(dir.Z) else inf
    local bx=if sx>0 then (cx+1)*cs else cx*cs
    local by=if sy>0 then (cy+1)*cs else cy*cs
    local bz=if sz>0 then (cz+1)*cs else cz*cs
    local tx=if sx~=0 then (bx-origin.X)/dir.X else inf
    local ty=if sy~=0 then (by-origin.Y)/dir.Y else inf
    local tz=if sz~=0 then (bz-origin.Z)/dir.Z else inf
    if tx<0 then tx=0 end;if ty<0 then ty=0 end;if tz<0 then tz=0 end
    local tcell=0
    while tcell<=bestT do
        -- Small bodies are stored by center cell while their AABB may touch an adjacent cell.
        -- Probe the 3x3x3 neighborhood and de-duplicate by body index.
        for ox=-1,1 do for oy=-1,1 do for oz=-1,1 do
            local key=packCell(cx+ox,cy+oy,cz+oz)
            if key then
                local i=heads[key] or 0
                while i~=0 do
                    if marks[i]~=epoch then
                        marks[i]=epoch
                        local t,n,p=self:_raycastBody(i,origin,dir,bestT,mask)
                        if t and t<bestT then bestT,bestI,bestN,bestP=t,i,n,p end
                    end
                    i=nextIn[i] or 0
                end
            end
        end end end
        local nextT=math.min(tx,ty,tz)
        if nextT>bestT or nextT==inf then break end
        local edge=nextT+1e-12
        if tx<=edge then cx+=sx;tx+=tdx end
        if ty<=edge then cy+=sy;ty+=tdy end
        if tz<=edge then cz+=sz;tz+=tdz end
        tcell=nextT
    end
    return bestI,bestT,bestP,bestN
end

function World:_sphereCastRaw(startPos,endPos,radius,mask,ignoreHandle)
    assert(typeof(startPos)=="vector" and typeof(endPos)=="vector" and type(radius)=="number" and radius>0,"physik.sphereCast: expected start/end vectors and radius > 0")
    mask=if mask==nil then 0xffffffff else mask
    local ignore=nil
    if ignoreHandle~=nil then ignore=self:_check(ignoreHandle) end
    local rvec=vector.create(radius,radius,radius)
    local mn=vector.create(math.min(startPos.X,endPos.X),math.min(startPos.Y,endPos.Y),math.min(startPos.Z,endPos.Z))-rvec
    local mx=vector.create(math.max(startPos.X,endPos.X),math.max(startPos.Y,endPos.Y),math.max(startPos.Z,endPos.Z))+rvec
    local handles,count=self:queryAabb(mn,mx,mask,self.castCandidates)
    self.castCandidates=handles
    local bestT,bestN,bestI=1,nil,nil
    for k=1,count do
        local i=self:_check(handles[k])
        if i==ignore then continue end
        local target=self.shape[i]
        local t,n
        if target.kind==SHAPE_SPHERE then
            t,n=sweepSphereSphere(startPos,endPos,radius,self.pos[i],target.radius)
        elseif target.kind==SHAPE_BOX then
            t,n=sweepSphereBox(startPos,endPos,radius,self.pos[i],self.qv[i],self.qw[i],target.halfExtents)
        elseif target.kind==SHAPE_CAPSULE then
            local a,b=capsuleEndpoints(self.pos[i],self.qv[i],self.qw[i],target.halfHeight)
            t,n=sweepSphereSegment(startPos,endPos,radius+target.radius,a,b)
        elseif target.kind==SHAPE_MESH then
            t,n=self:_sweepSphereMesh(i,startPos,endPos,radius)
        elseif target.kind==SHAPE_HEIGHTFIELD then
            t,n=self:_sweepSphereHeightfield(i,startPos,endPos,radius)
        elseif target.kind==SHAPE_COMPOUND then
            t,n=self:_sweepSphereCompound(i,startPos,endPos,radius)
        end
        if t and t<bestT then bestT,bestN,bestI=t,n,i end
    end
    if not bestI then return nil end
    local center=startPos+(endPos-startPos)*bestT
    return bestI,bestT,bestN,center
end

function World:sphereCast(startPos,endPos,radius,mask,ignoreHandle)
    local i,t,n,center=self:_sphereCastRaw(startPos,endPos,radius,mask,ignoreHandle)
    if not i then return nil end
    local motion=endPos-startPos
    local distance=vector.magnitude(motion)*t
    return {handle=makeHandle(i,self.generation[i]),fraction=t,distance=distance,center=center,point=center-n*radius,normal=n,userData=self.userData[i]}
end

function World:sphereCastBatchSoA(starts,ends,radii,masks,ignoreHandles,out)
    assert(type(starts)=="table" and type(ends)=="table" and #starts==#ends,"physik.sphereCastBatchSoA: starts and ends must be equal-length arrays")
    out=out or {}
    local hit=out.hit or {};out.hit=hit
    local handles=out.handles or {};out.handles=handles
    local fractions=out.fractions or {};out.fractions=fractions
    local distances=out.distances or {};out.distances=distances
    local centers=out.centers or {};out.centers=centers
    local points=out.points or {};out.points=points
    local normals=out.normals or {};out.normals=normals
    local userData=out.userData or {};out.userData=userData
    local count=#starts
    local oldCount=out.count or 0
    for q=1,count do
        local radius=if type(radii)=="table" then radii[q] else radii
        local mask=if masks then masks[q] else nil
        local ignore=if ignoreHandles then ignoreHandles[q] else nil
        local i,t,n,center=self:_sphereCastRaw(starts[q],ends[q],radius,mask,ignore)
        if i then
            local motion=ends[q]-starts[q]
            hit[q]=true;handles[q]=makeHandle(i,self.generation[i]);fractions[q]=t;distances[q]=vector.magnitude(motion)*t
            centers[q]=center;points[q]=center-n*radius;normals[q]=n;userData[q]=self.userData[i]
        else
            hit[q]=false;handles[q]=nil;fractions[q]=nil;distances[q]=nil;centers[q]=nil;points[q]=nil;normals[q]=nil;userData[q]=nil
        end
    end
    for q=count+1,oldCount do
        hit[q]=nil;handles[q]=nil;fractions[q]=nil;distances[q]=nil;centers[q]=nil;points[q]=nil;normals[q]=nil;userData[q]=nil
    end
    out.count=count
    return out,count
end

function World:raycast(origin,direction,maxDistance,mask)
    local i,t,p,n=self:_raycastRaw(origin,direction,maxDistance,mask)
    if not i then return nil end
    return { handle=makeHandle(i,self.generation[i]), distance=t, point=p, normal=n, userData=self.userData[i] }
end

function World:raycastBatch(origins,directions,maxDistances,masks,out)
    assert(type(origins)=="table" and type(directions)=="table" and #origins==#directions,"physik.raycastBatch: origins and directions must be equal-length arrays")
    out=out or {}
    local count=#origins
    local oldCount=out.count or #out
    for q=1,count do
        local i,t,p,n=self:_raycastRaw(origins[q],directions[q],if maxDistances then maxDistances[q] else nil,if masks then masks[q] else nil)
        if i then
            local hit=out[q]
            if type(hit)~="table" then hit={};out[q]=hit end
            hit.handle=makeHandle(i,self.generation[i])
            hit.distance=t
            hit.point=p
            hit.normal=n
            hit.userData=self.userData[i]
        else
            out[q]=nil
        end
    end
    for q=count+1,oldCount do out[q]=nil end
    out.count=count
    return out
end

function World:raycastBatchSoA(origins,directions,maxDistances,masks,out)
    assert(type(origins)=="table" and type(directions)=="table" and #origins==#directions,"physik.raycastBatchSoA: origins and directions must be equal-length arrays")
    out=out or {}
    local hit=out.hit or {};out.hit=hit
    local handles=out.handles or {};out.handles=handles
    local distances=out.distances or {};out.distances=distances
    local points=out.points or {};out.points=points
    local normals=out.normals or {};out.normals=normals
    local userData=out.userData or {};out.userData=userData
    local count=#origins
    local oldCount=out.count or 0
    for q=1,count do
        local i,t,p,n=self:_raycastRaw(origins[q],directions[q],if maxDistances then maxDistances[q] else nil,if masks then masks[q] else nil)
        if i then
            hit[q]=true
            handles[q]=makeHandle(i,self.generation[i])
            distances[q]=t
            points[q]=p
            normals[q]=n
            userData[q]=self.userData[i]
        else
            hit[q]=false
            handles[q]=nil
            distances[q]=nil
            points[q]=nil
            normals[q]=nil
            userData[q]=nil
        end
    end
    for q=count+1,oldCount do
        hit[q]=nil;handles[q]=nil;distances[q]=nil;points[q]=nil;normals[q]=nil;userData[q]=nil
    end
    out.count=count
    return out,count
end


local function bodyOverlapsAabb(self,i,mn,mx)
    local shape=self.shape[i]
    if shape.kind==SHAPE_SPHERE then
        local p=self.pos[i];local r=shape.radius+self.contactMargin
        return not (p.X+r<mn.X or p.X-r>mx.X or p.Y+r<mn.Y or p.Y-r>mx.Y or p.Z+r<mn.Z or p.Z-r>mx.Z)
    end
    local a,b=self.aabbMin[i],self.aabbMax[i]
    return not (b.X<mn.X or a.X>mx.X or b.Y<mn.Y or a.Y>mx.Y or b.Z<mn.Z or a.Z>mx.Z)
end

function World:_queryAabbScan(mn,mx,mask,out)
    local n=0
    for i=1,self.high do
        if not self.alive[i] or bit32.band(mask,self.group[i])==0 then continue end
        if not bodyOverlapsAabb(self,i,mn,mx) then continue end
        n+=1
        out[n]=makeHandle(i,self.generation[i])
    end
    return n
end

function World:queryAabb(mn,mx,mask,out)
    assert(typeof(mn)=="vector" and typeof(mx)=="vector","physik.queryAabb: bounds must be vectors")
    mask=if mask==nil then 0xffffffff else mask
    out=out or {}
    self:_ensureQueryGrid()
    local cs=self.broadphaseCellSize
    -- Small bodies live in their center cell, but may extend into the immediately adjacent cell.
    local x0,y0,z0=math.floor(mn.X/cs)-1,math.floor(mn.Y/cs)-1,math.floor(mn.Z/cs)-1
    local x1,y1,z1=math.floor(mx.X/cs)+1,math.floor(mx.Y/cs)+1,math.floor(mx.Z/cs)+1
    local cells=(x1-x0+1)*(y1-y0+1)*(z1-z0+1)
    local n=0
    if cells>self.queryAabbCellLimit then
        n=self:_queryAabbScan(mn,mx,mask,out)
    else
        -- Oversized bodies are outside the grid and are normally few.
        for k=1,#self.gridLarge do
            local i=self.gridLarge[k]
            if self.alive[i] and bit32.band(mask,self.group[i])~=0 then
                if bodyOverlapsAabb(self,i,mn,mx) then
                    n+=1;out[n]=makeHandle(i,self.generation[i])
                end
            end
        end
        self.queryMarkEpoch+=1
        local epoch=self.queryMarkEpoch
        local marks,heads,nextIn=self.queryMarks,self.gridHeads,self.gridNext
        for x=x0,x1 do for y=y0,y1 do for z=z0,z1 do
            local key=packCell(x,y,z)
            if key then
                local i=heads[key] or 0
                while i~=0 do
                    if marks[i]~=epoch then
                        marks[i]=epoch
                        if bit32.band(mask,self.group[i])~=0 then
                            if bodyOverlapsAabb(self,i,mn,mx) then
                                n+=1;out[n]=makeHandle(i,self.generation[i])
                            end
                        end
                    end
                    i=nextIn[i] or 0
                end
            end
        end end end
    end
    for i=#out,n+1,-1 do out[i]=nil end
    return out,n
end


-- Ephemeral zero-copy view of hot linear state for trusted bulk consumers. The returned arrays
-- are Physik-owned and must be treated as read-only; their contents may change on the next world
-- mutation/step. This avoids copying positions/velocities into an intermediate drain when a caller
-- is immediately streaming them into another packed store.
function World:linearStateView()
    local view=self.linearView
    view.high=self.high
    view.dynamicCode=BODY_DYNAMIC
    return view
end

function World:writeLinearState(out,opts)
    out=out or {}
    opts=opts or {}
    local includeHandles=opts.handles~=false
    local includeUserData=opts.userData==true
    local handles=if includeHandles then (out.handles or {}) else nil
    local userData=if includeUserData then (out.userData or {}) else nil
    local positions=out.positions or {}
    local velocities=out.velocities or {}
    out.handles=handles;out.userData=userData;out.positions=positions;out.velocities=velocities
    local n=0
    local awakeOnly=opts.awakeOnly==true
    local dynamicOnly=opts.dynamicOnly==true
    for i=1,self.high do
        if not self.alive[i] then continue end
        if awakeOnly and self.sleeping[i] then continue end
        if dynamicOnly and self.bodyType[i]~=BODY_DYNAMIC then continue end
        n+=1
        if handles then handles[n]=makeHandle(i,self.generation[i]) end
        if userData then userData[n]=self.userData[i] end
        positions[n]=self.pos[i]
        velocities[n]=self.linVel[i]
    end
    if handles then for i=#handles,n+1,-1 do handles[i]=nil end end
    if userData then for i=#userData,n+1,-1 do userData[i]=nil end end
    for i=#positions,n+1,-1 do positions[i]=nil;velocities[i]=nil end
    out.count=n
    return out,n
end

function World:writeTransforms(out,opts)
    out=out or {}
    opts=opts or {}
    local handles=out.handles or {};out.handles=handles
    local positions=out.positions or {};out.positions=positions
    local rotations=out.rotations or {};out.rotations=rotations
    local rotationW=out.rotationW or {};out.rotationW=rotationW
    local velocities=out.velocities or {};out.velocities=velocities
    local angularVelocities=out.angularVelocities or {};out.angularVelocities=angularVelocities
    local n=0
    local awakeOnly=opts.awakeOnly==true
    local dynamicOnly=opts.dynamicOnly==true
    for i=1,self.high do
        if not self.alive[i] then continue end
        if awakeOnly and self.sleeping[i] then continue end
        if dynamicOnly and self.bodyType[i]~=BODY_DYNAMIC then continue end
        n+=1
        handles[n]=makeHandle(i,self.generation[i])
        positions[n]=self.pos[i]
        rotations[n]=self.qv[i]
        rotationW[n]=self.qw[i]
        velocities[n]=self.linVel[i]
        angularVelocities[n]=self.angVel[i]
    end
    for i=#handles,n+1,-1 do handles[i]=nil;positions[i]=nil;rotations[i]=nil;rotationW[i]=nil;velocities[i]=nil;angularVelocities[i]=nil end
    out.count=n
    return out,n
end


function World:applyTransforms(data,wake)
    assert(type(data)=="table" and type(data.handles)=="table" and type(data.positions)=="table","physik.applyTransforms: expected writeTransforms-compatible data")
    local count=data.count or #data.handles
    local qvs,qws=data.rotations,data.rotationW
    local vels,ang=data.velocities,data.angularVelocities
    for n=1,count do
        local i=self:_check(data.handles[n])
        self.pos[i]=data.positions[n]
        if qvs and qws and qvs[n] and qws[n] then
            self.qv[i],self.qw[i]=qNormalize(qvs[n],qws[n])
            self.axis0[i],self.axis1[i],self.axis2[i]=qAxes(self.qv[i],self.qw[i])
        end
        if vels and vels[n] then self.linVel[i]=vels[n] end
        if ang and ang[n] then self.angVel[i]=ang[n] end
        self.prevPos[i]=self.pos[i]
        self:_updateAabb(i)
        if wake~=false and self.bodyType[i]==BODY_DYNAMIC then self.sleeping[i]=false;self.sleepTime[i]=0 end
    end
    -- Topology is unchanged, but authoritative teleports invalidate persistent anchors and query cells.
    self.cacheEpoch+=2
    self.gridValid=false
    return count
end

function World:writeContacts(out)
    out=out or {}
    local aa=out.a or {};out.a=aa
    local bb=out.b or {};out.b=bb
    local normals=out.normals or {};out.normals=normals
    local points=out.points or {};out.points=points
    local penetration=out.penetration or {};out.penetration=penetration
    local impulses=out.normalImpulse or {};out.normalImpulse=impulses
    local sensors=out.sensor or {};out.sensor=sensors
    local n=self.contactCount
    for c=1,n do
        local a,b=self.ca[c],self.cb[c]
        aa[c]=makeHandle(a,self.generation[a]);bb[c]=makeHandle(b,self.generation[b])
        normals[c]=self.cn[c];points[c]=self.cp[c];penetration[c]=self.cpen[c];impulses[c]=self.cnormalImpulse[c]
        sensors[c]=self.sensor[a] or self.sensor[b]
    end
    for c=#aa,n+1,-1 do aa[c]=nil;bb[c]=nil;normals[c]=nil;points[c]=nil;penetration[c]=nil;impulses[c]=nil;sensors[c]=nil end
    out.count=n
    return out,n
end

function World:writeContactEvents(out,includeStay)
    out=out or {}
    local aa=out.a or {};out.a=aa
    local bb=out.b or {};out.b=bb
    local kind=out.kind or {};out.kind=kind
    local sensors=out.sensor or {};out.sensor=sensors
    local n=0
    for e=1,self.contactEventCount do
        local k=self.eventKind[e]
        if includeStay~=false or k~=EVENT_STAY then
            local a,b=self.eventA[e],self.eventB[e]
            if self.alive[a] and self.alive[b] then
                n+=1
                aa[n]=makeHandle(a,self.generation[a]);bb[n]=makeHandle(b,self.generation[b])
                kind[n]=k;sensors[n]=self.eventSensor[e]
            end
        end
    end
    for i=#aa,n+1,-1 do aa[i]=nil;bb[i]=nil;kind[i]=nil;sensors[i]=nil end
    out.count=n
    return out,n
end

function World:advance(dt,fixedDt,maxSteps)
    assert(type(dt)=="number" and dt>=0,"physik.advance: dt must be >= 0")
    local h=fixedDt or self.fixedDt
    local limit=maxSteps or self.maxCatchUpSteps
    assert(h>0 and limit>=1 and limit%1==0,"physik.advance: fixedDt must be > 0 and maxSteps a positive integer")
    self.accumulator+=dt
    local wanted=math.floor(self.accumulator/h)
    local steps=math.min(wanted,limit)
    for _=1,steps do self:step(h) end
    self.accumulator-=steps*h
    local dropped=false
    if self.accumulator>=h then self.accumulator=self.accumulator%h;dropped=true end
    return steps,self.accumulator/h,dropped
end

function World:contacts()
    local out=table.create(self.contactCount)
    for c=1,self.contactCount do
        local a,b=self.ca[c],self.cb[c]
        out[c]={
            a=makeHandle(a,self.generation[a]),
            b=makeHandle(b,self.generation[b]),
            normal=self.cn[c],point=self.cp[c],penetration=self.cpen[c],
            normalImpulse=self.cnormalImpulse[c],sensor=self.sensor[a] or self.sensor[b],
        }
    end
    return out
end

function World:stats()
    local asleep,sensorContacts=0,0
    for i=1,self.high do if self.alive[i] and self.sleeping[i] then asleep+=1 end end
    for c=1,self.contactCount do if self.sensor[self.ca[c]] or self.sensor[self.cb[c]] then sensorContacts+=1 end end
    return {
        frame=self.frame,bodies=self.aliveCount,sleeping=asleep,islands=self.islandCount,
        candidates=self.candidateCount,contacts=self.contactCount,sensorContacts=sensorContacts,activeContacts=self.activeContactCount,manifolds=self.manifoldCount,warmStarted=self.warmStartedCount,ccdHits=self.ccdHitCount,
        joints=self.jointCount,
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
    w.trackContactEvents=opts.contactEvents==true
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
    w.meshQueryCellLimit=opts.meshQueryCellLimit or 4096
    assert(w.meshQueryCellLimit>=1 and w.meshQueryCellLimit%1==0,"physik.newWorld: meshQueryCellLimit must be a positive integer")
    w.queryRayCellLimit=opts.queryRayCellLimit or 2048
    assert(w.queryRayCellLimit>=1 and w.queryRayCellLimit%1==0,"physik.newWorld: queryRayCellLimit must be a positive integer")
    w.queryAabbCellLimit=opts.queryAabbCellLimit or 8192
    assert(w.queryAabbCellLimit>=1 and w.queryAabbCellLimit%1==0,"physik.newWorld: queryAabbCellLimit must be a positive integer")
    w.fixedDt=opts.fixedDt or (1/60)
    w.maxCatchUpSteps=opts.maxCatchUpSteps or 8
    w.accumulator=0
    w.gridHeads={};w.gridKeys={};w.gridNext={};w.gridLarge={};w.gridX={};w.gridY={};w.gridZ={};w.gridKey={}
    w.gridKeyX={};w.gridKeyY={};w.gridKeyZ={}
    w.gridOccupied = buffer.create(GridFilterSize)
    w.gridDenseKeys={};w.gridPairA={};w.gridPairB={};w.gridDenseCount=0;w.gridPairCount=0
    w.gridTopologyReady=false;w.gridTopologyDirty=true
    w.gridValid=false;w.gridSmallCount=0;w.queryMarks={};w.queryMarkEpoch=0;w.castCandidates={}
    w.frame=0
    w.cacheEpoch=0
    w.ccdHitCount=0
    w.ccdCount=0
    w.high=0
    w.aliveCount=0
    w.free={}
    w.order={}
    w.generation={}
    w.alive={}
    w.bodyType={}
    w.shape={}
    w.sphereRadius = {}
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
    w.sensor={}
    w.ccd={};w.prevPos={}
    w.linearView={alive=w.alive,bodyType=w.bodyType,userData=w.userData,positions=w.pos,velocities=w.linVel,high=0,dynamicCode=BODY_DYNAMIC}
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
    w.contactEventEpoch=0;w.contactEventCount=0
    w.contactPairSeen={};w.contactPairActive={};w.contactPairA={};w.contactPairB={};w.contactPairSensor={}
    w.contactPairPrev={};w.contactPairNext={}
    w.eventA={};w.eventB={};w.eventKind={};w.eventSensor={}
    w.clipA={};w.clipB={};w.manifoldPoints={};w.manifoldDepths={}
    w.geomNormal={};w.geomPoint={};w.geomPen={};w.geomContactCount=0
    w.meshMarks={};w.meshMarkEpoch=0
    w.islandParent={};w.islandSize={};w.islandAwake={};w.islandRoots={};w.islandCanSleep={};w.islandMinTime={}
    w.joints={};w.jointNoCollide={};w.nextJoint=1;w.jointCount=0
    return w
end

M.contactEvent = table.freeze({ enter=EVENT_ENTER, stay=EVENT_STAY, exit=EVENT_EXIT })
M._constants = table.freeze({ sphere=SHAPE_SPHERE, box=SHAPE_BOX, capsule=SHAPE_CAPSULE, mesh=SHAPE_MESH, heightfield=SHAPE_HEIGHTFIELD, compound=SHAPE_COMPOUND })
return table.freeze(M)
