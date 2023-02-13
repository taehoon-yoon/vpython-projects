import vpython as vp
import numpy as np


def get_idx(num, depth, row, col):
    return num * num * depth + num * row + col


def vpVec2npVec(vec):
    return np.array([vec.x, vec.y, vec.z])


def cal_acc(pos, ci, cj, spring_k, spring_L0, mass_mat, gravity_acc, velocity, drag_coef):
    """
    Calculating acceleration on each particle(node)
    1. hook force (spring force)
    2. gravity
    3. drag force
    """
    acc = np.zeros_like(pos)
    # 1. hook force calculation
    displacement_vec = pos[ci, :] - pos[cj, :]
    displacement_len = np.linalg.norm(displacement_vec, axis=1)
    delta_len = displacement_len - spring_L0
    hook_force = -spring_k[..., None] * delta_len[..., None] * displacement_vec / displacement_len[..., None]
    np.add.at(acc, ci, hook_force)
    np.add.at(acc, cj, -hook_force)
    # 2. gravity
    acc[:, 1] += -gravity_acc * mass_mat
    # 3. drag force
    acc -= drag_coef * velocity
    return acc


def get_box_min_max(box_shape):
    """
    Get box lower left corner and upper right corner in object coordinate(box coordinate)
    Assuming box's center is located at origin
    """
    min_ = np.array([-box_shape.x / 2, -box_shape.y / 2, -box_shape.z / 2])
    max_ = np.array([box_shape.x / 2, box_shape.y / 2, box_shape.z / 2])
    return min_, max_


def pointInAABB(min_, max_, point):
    """
    AABB algorithm to check if the point(node) is inside the box if point is inside the box
    then collision has occured.
    """
    # point is numpy array (# of node, 3)
    is_inside = np.stack([(min_ < point).all(axis=1), (max_ > point).all(axis=1)], axis=0).all(axis=0)
    return is_inside


def get_camera_extrinsic(camera_x, camera_y, camera_origin):
    """
    Given camera's origin and camera's x and y axis, calculate camera extrinsic matrix.
    w2c is world to camera transformation matrix.
    c2w is camera to world transformation matrix.
    camera_x and camera_y and camera_origin are all in world coordinate sysytem.
    w2c, c2w all 4X4 matrix since we use homogeneous coordinate
    """
    # type conversion
    camera_origin = vpVec2npVec(camera_origin)
    camera_x = vpVec2npVec(camera_x)
    camera_y = vpVec2npVec(camera_y)

    # Normalize
    camera_x /= np.linalg.norm(camera_x)
    camera_y /= np.linalg.norm(camera_y)
    camera_z = np.cross(camera_x, camera_y)
    camera_z /= np.linalg.norm(camera_z)

    rotation_mat = np.vstack([camera_x, camera_y, camera_z])
    rotation_matT = rotation_mat.T
    translation = -np.matmul(rotation_mat, camera_origin)

    w2c = np.concatenate([rotation_mat, translation[..., None]], axis=1)
    w2c = np.concatenate([w2c, np.array([[0, 0, 0, 1]])])

    c2w = np.concatenate([rotation_matT, camera_origin[..., None]], axis=1)
    c2w = np.concatenate([c2w, np.array([[0, 0, 0, 1]])])
    return w2c, c2w


def collision(objects, pos, vel):
    """
    collision detection algorithm, currently only box object is supported.
    """
    for obj in objects:
        if obj.type == 'box':
            """
            Boxes are in general not aligned with world coordinate. If boxes are all aligned with
            world coordinate, what we have to do is just apply AABB algorithm. But since box can be tilted and rotated,
            we cannot simply apply AABB algorithm to such boxes to detect collision. 
            So in the case boxes are not aligned, we transform all the points (which are represented in world coord)
            to box coordinate frame using the box w2c matrix. Thus we get points in box coordinate. And since in
            box coordinate each boxes are all aligned, we can use AABB algorithm to detect collision.
            """
            w2c = obj.w2c
            c2w = obj.c2w
            homogeneous_coord = np.concatenate([pos, np.ones([pos.shape[0], 1])], axis=1).T  # (4, # of nodes)
            pos_camera = np.matmul(w2c, homogeneous_coord)[:3, :].T  # (# of node, 3)
            # Now in camera coord, box center is in the origin, axis points x direction, up points y direction
            is_inside = pointInAABB(obj.min, obj.max, pos_camera)  # (# of nodes, )

            """
            from is_inside, we can know what particles collide with box. For the collided particles 
            we flip the particles position based on the upper side of the box.
            """
            pos_camera[is_inside, 1] = obj.size.y - pos_camera[is_inside, 1]

            # Now back to world coord
            homogeneous_coord_back = np.concatenate([pos_camera, np.ones([pos_camera.shape[0], 1])], axis=1).T
            pos = np.matmul(c2w, homogeneous_coord_back)[:3, :].T # (# of node, 3)
            # Now particles are now in world coordinate with collided particles flipped based on upper side of box

            # flip velocity, also we have to flip velocity vector for those collided with box
            # we assume specular reflection, -2(v dot n)n + v
            normal_vec = vpVec2npVec(obj.up)
            normal_vec = normal_vec / np.linalg.norm(normal_vec)
            vel[is_inside] = (-2 * np.matmul(vel[is_inside], normal_vec))[:, None] * normal_vec[None, :] + vel[
                is_inside]
    return pos, vel


N = 5  # number of particles(nodes) along each axis
dt = 1 / 30
spring_coef = 150
gravitational_acceleration = 9.8
node_mass = 0.1
drag_coef = 0.3

vp.scene = vp.canvas(title='spring network',
                     width=1700, height=750,
                     center=vp.vector(2, 0, 0.1), background=vp.color.cyan)

# Making mass spring system (particles attached with each other with spring), we will add springs below
xlin = np.linspace(6, 8, N)
ylin = np.linspace(2, 4, N)
zlin = np.linspace(-1, 1, N)
x, y, z = np.meshgrid(xlin, ylin, zlin)
x = x.flatten()
y = y.flatten()
z = z.flatten()
pos = np.vstack([x, y, z]).T
acc = np.zeros_like(pos)
velocity = 0.01 * np.random.randn(x.shape[0], 3)  # random initial velocity for each particles.
node_array = []  # node == particle
for node in pos:
    node_ = vp.sphere(pos=vp.vector(*node), radius=0.1, color=vp.color.blue)
    node_.mass = node_mass
    node_array.append(node_)

# adding objects(which can collide with particles) including ground etc
object_array = []
ground1 = vp.box(pos=vp.vec(4, -1, 0), size=vp.vec(3, 0.2, 6), color=vp.color.green)
ground1.type = 'box'
ground1.w2c, ground1.c2w = get_camera_extrinsic(ground1.axis, ground1.up, ground1.pos)  # 4X4
ground1.min, ground1.max = get_box_min_max(ground1.size)
object_array.append(ground1)

clamp = vp.box(pos=vp.vec(4, -1, 0), size=vp.vec(10, 0.2, 6), color=vp.color.green)
clamp.type = 'box'
clamp.up = vp.vec(-1, 2, 0)
clamp.w2c, clamp.c2w = get_camera_extrinsic(clamp.axis, clamp.up, clamp.pos)  # 4X4
clamp.min, clamp.max = get_box_min_max(clamp.size)
object_array.append(clamp)

ground2 = vp.box(pos=vp.vec(2, -6, 0), size=vp.vec(20, 0.2, 6), color=vp.color.green)
ground2.type = 'box'
ground2.w2c, ground2.c2w = get_camera_extrinsic(ground2.axis, ground2.up, ground2.pos)  # 4X4
ground2.min, ground2.max = get_box_min_max(ground2.size)
object_array.append(ground2)

# Adding springs
ci = []
cj = []

"""
O--O    (O is particle -- is spring)
"""
for d in range(N):
    for r in range(N):
        for c in range(N - 1):
            ci.append(get_idx(N, d, r, c))
            cj.append(get_idx(N, d, r, c + 1))

"""
O
|
|
O
"""
for d in range(N):
    for r in range(N - 1):
        for c in range(N):
            ci.append(get_idx(N, d, r, c))
            cj.append(get_idx(N, d, r + 1, c))

"""
O
 \
  \
   O
"""
for d in range(N):
    for r in range(N - 1):
        for c in range(N - 1):
            ci.append(get_idx(N, d, r, c))
            cj.append(get_idx(N, d, r + 1, c + 1))

"""
   O
  /
 /
O
"""
for d in range(N):
    for r in range(N - 1):
        for c in range(N - 1):
            ci.append(get_idx(N, d, r + 1, c))
            cj.append(get_idx(N, d, r, c + 1))


"""
Doing exactly same thing like above 4 for loop. This time connecting in depth dimension.
"""
for d in range(N - 1):
    for r in range(N):
        for c in range(N):
            ci.append(get_idx(N, d, r, c))
            cj.append(get_idx(N, d + 1, r, c))

for d in range(N - 1):
    for r in range(N):
        for c in range(N - 1):
            ci.append(get_idx(N, d, r, c))
            cj.append(get_idx(N, d + 1, r, c + 1))

for d in range(N - 1):
    for r in range(N):
        for c in range(N - 1):
            ci.append(get_idx(N, d, r, c + 1))
            cj.append(get_idx(N, d + 1, r, c))

for d in range(N - 1):
    for r in range(N-1):
        for c in range(N):
            ci.append(get_idx(N, d, r, c))
            cj.append(get_idx(N, d + 1, r+1, c))

for d in range(N - 1):
    for r in range(N-1):
        for c in range(N):
            ci.append(get_idx(N, d, r+1, c))
            cj.append(get_idx(N, d + 1, r, c))

edge_array = []
for i, j in zip(ci, cj):
    spring = vp.cylinder(pos=node_array[i].pos, axis=node_array[j].pos - node_array[i].pos, radius=0.01)
    spring.k = spring_coef
    edge_array.append(spring)

spring_L0 = np.array([edge.axis.mag for edge in edge_array])
spring_k = np.array([edge.k for edge in edge_array])
mass_mat = np.array([node.mass for node in node_array])

while True:
    vp.rate(1 / dt)

    # Leapfrog integration
    velocity += acc * dt / 2
    pos += velocity * dt
    pos, velocity = collision(object_array, pos, velocity)
    acc = cal_acc(pos, ci, cj, spring_k, spring_L0, mass_mat, gravitational_acceleration, velocity, drag_coef)
    velocity += acc * dt / 2

    # update particle position for rendering
    for node, updated_pos in zip(node_array, pos):
        node.pos = vp.vec(*updated_pos)

    # update spring position for rendering
    for i, j, spring in zip(ci, cj, edge_array):
        spring.pos = node_array[i].pos
        spring.axis = node_array[j].pos - node_array[i].pos
