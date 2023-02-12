import vpython as vp
import numpy as np


def get_idx(num, depth, row, col):
    return num * num * depth + num * row + col


def vpVec2npVec(vec):
    return np.array([vec.x, vec.y, vec.z])


def cal_acc(pos, ci, cj, spring_k, spring_L0, mass_mat, gravity_acc, velocity, drag_coef):
    acc = np.zeros_like(pos)
    displacement_vec = pos[ci, :] - pos[cj, :]
    displacement_len = np.linalg.norm(displacement_vec, axis=1)
    delta_len = displacement_len - spring_L0
    hook_force = -spring_k[..., None] * delta_len[..., None] * displacement_vec / displacement_len[..., None]
    np.add.at(acc, ci, hook_force)
    np.add.at(acc, cj, -hook_force)
    acc[:, 1] += -gravity_acc * mass_mat
    acc -= drag_coef * velocity
    return acc


def get_box_min_max(box_shape):
    min_ = np.array([-box_shape.x / 2, -box_shape.y / 2, -box_shape.z / 2])
    max_ = np.array([box_shape.x / 2, box_shape.y / 2, box_shape.z / 2])
    return min_, max_


def pointInAABB(min_, max_, point):
    # point is numpy array (# of node, 3)
    is_inside = np.stack([(min_ < point).all(axis=1), (max_ > point).all(axis=1)], axis=0).all(axis=0)
    return is_inside


def get_camera_extrinsic(camera_x, camera_y, camera_origin):
    camera_origin = vpVec2npVec(camera_origin)
    camera_x = vpVec2npVec(camera_x)
    camera_y = vpVec2npVec(camera_y)
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
    for obj in objects:
        if obj.type == 'box':
            w2c = obj.w2c
            c2w = obj.c2w
            homogeneous_coord = np.concatenate([pos, np.ones([pos.shape[0], 1])], axis=1).T
            pos_camera = np.matmul(w2c, homogeneous_coord)[:3, :].T  # (# of node, 3)
            # Now in camera coord, box pos is the origin, axis points x direction, up points y direction
            is_inside = pointInAABB(obj.min, obj.max, pos_camera)  # (# of nodes, )
            pos_camera[is_inside, 1] = obj.size.y - pos_camera[is_inside, 1]  # flip node
            # Now back to world coord
            homogeneous_coord_back = np.concatenate([pos_camera, np.ones([pos_camera.shape[0], 1])], axis=1).T
            pos = np.matmul(c2w, homogeneous_coord_back)[:3, :].T

            # flip velocity
            # -2(v dot n)n + v
            normal_vec = vpVec2npVec(obj.up)
            normal_vec = normal_vec / np.linalg.norm(normal_vec)
            vel[is_inside] = (-2 * np.matmul(vel[is_inside], normal_vec))[:, None] * normal_vec[None, :] + vel[
                is_inside]
    return pos, vel


N = 5
dt = 1 / 30
spring_coef = 150
gravitational_acceleration = 9.8
node_mass = 0.1
drag_coef = 0.3

vp.scene = vp.canvas(title='spring network',
                     width=1700, height=750,
                     center=vp.vector(2, 0, 0.1), background=vp.color.cyan)
xlin = np.linspace(6, 8, N)
ylin = np.linspace(2, 4, N)
zlin = np.linspace(-1, 1, N)
x, y, z = np.meshgrid(xlin, ylin, zlin)
x = x.flatten()
y = y.flatten()
z = z.flatten()
pos = np.vstack([x, y, z]).T
acc = np.zeros_like(pos)
velocity = 0.01 * np.random.randn(x.shape[0], 3)
node_array = []
for node in pos:
    node_ = vp.sphere(pos=vp.vector(*node), radius=0.1, color=vp.color.blue)
    node_.mass = node_mass
    node_array.append(node_)

object_array = []
ground1 = vp.box(pos=vp.vec(4, -1, 0), size=vp.vec(3, 0.2, 6), color=vp.color.green)
ground1.type = 'box'
# ground1.up=vp.vec(-1,0.5,0)
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

ci = []
cj = []

for d in range(N):
    for r in range(N):
        for c in range(N - 1):
            ci.append(get_idx(N, d, r, c))
            cj.append(get_idx(N, d, r, c + 1))

for d in range(N):
    for r in range(N - 1):
        for c in range(N):
            ci.append(get_idx(N, d, r, c))
            cj.append(get_idx(N, d, r + 1, c))

for d in range(N):
    for r in range(N - 1):
        for c in range(N - 1):
            ci.append(get_idx(N, d, r, c))
            cj.append(get_idx(N, d, r + 1, c + 1))
for d in range(N):
    for r in range(N - 1):
        for c in range(N - 1):
            ci.append(get_idx(N, d, r + 1, c))
            cj.append(get_idx(N, d, r, c + 1))

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

    velocity += acc * dt / 2
    pos += velocity * dt
    pos, velocity = collision(object_array, pos, velocity)
    acc = cal_acc(pos, ci, cj, spring_k, spring_L0, mass_mat, gravitational_acceleration, velocity, drag_coef)
    velocity += acc * dt / 2

    for node, updated_pos in zip(node_array, pos):
        node.pos = vp.vec(*updated_pos)

    for i, j, spring in zip(ci, cj, edge_array):
        spring.pos = node_array[i].pos
        spring.axis = node_array[j].pos - node_array[i].pos
