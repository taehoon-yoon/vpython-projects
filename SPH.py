import vpython as vp
import numpy as np
from scipy.special import gamma
import matplotlib
import matplotlib.cm as cm


def kernel(pos_arr, h):
    """
    :param pos_arr: shape=[..., 3] / in particular [n, m, 3]
    :param h:
    :return: gaussian smoothing kernel value / shape=[...] / in particular [n, m]
    """
    distance = np.sum(pos_arr ** 2, axis=-1)
    return 1 / (h ** 3 * np.pi ** (3 / 2)) * np.exp(-distance / (h ** 2))


def kernel_derivatice(pos_arr, h):
    """
    :param pos_arr: shape=[...,3] / in particular [n, m, 3]
    :param h:
    :return: derivative of gaussign smoothing kernel at each position / shape=[..., 3] / in particular [n, m, 3]
    """
    distance = np.sum(pos_arr ** 2, axis=-1, keepdims=True)
    return -2 / (h ** 5 * np.pi ** (3 / 2)) * np.exp(-distance / (h ** 2)) * pos_arr


def get_pairwise_seperation(ri, rj):
    """
    Get pairwise seperation between 2 set of coordinate vector (ri - rj)
    :param ri: shape=[n, 3]
    :param rj: shape=[m, 3]
    :return: shape=[n, m, 3]
    """
    return ri[:, None, :] - rj[None, ...]


def cal_density(samping_loc, SPH_particle_loc, m, h):
    """
    Get Density at sampling loctions from SPH particle distribution
    :param samping_loc: samping location where we ant to estimate density using SPH formula / shape=[n,3]
    :param SPH_particle_loc: particle location / shape=[m,3]
    :param m: mass
    :param h: smoothing length
    :return: density at each sampled location / shape=[n,1]
    """
    pairwise_sep = get_pairwise_seperation(samping_loc, SPH_particle_loc)
    density = np.sum(m * kernel(pairwise_sep, h), axis=-1, keepdims=True)  # [n, 1]
    return density


def cal_pressure(density, k, n):
    """
    Equation of state  P=k * rho^(1+1/n)  -> polytropic equation of state
    :param density: density rho / shape=[n, 1]
    :param k: eqation of state constant
    :param n: polytropic index
    :return: [n,1]
    """
    return k * (density ** (1.0 + 1.0 / n))


def cal_acc(pos, vel, h, Lambda, nu, m, k, n):
    """
    calculate acceleration of each particle using SPH approximation
    :param pos: [N, 3]
    :param vel: [N, 3]
    :param h: smoothing length
    :param Lambda: external force constant (gravitational force)
    :param nu: viscosity for damping motion
    :param m: mass
    :param k: equation of state constant
    :param n: polytropic constant
    :return: [N, 3]
    """
    density = cal_density(pos, pos, m, h)  # [N,1]
    pressure = cal_pressure(density, k, n)  # [N,1]
    grad_kernel = kernel_derivatice(get_pairwise_seperation(pos, pos), h)  # [N,N,3]
    acc = -np.sum(m * ((pressure / (density ** 2)) + (pressure / (density ** 2)).T)[..., None] * grad_kernel,
                  axis=1)  # [N,3]
    acc += - Lambda * pos - nu * vel
    return acc


def main():
    vp.scene = vp.canvas(title='SPH Star Simulation', width=750, height=750)
    dt = 1 / 60
    N = 400
    M = 2
    m = M / N
    R = 0.75
    r = 0.02
    h = 0.1
    k = 0.1
    n = 1
    nu = 1
    Lambda = 2 * k * (1 + n) * np.pi ** (-3 / (2 * n)) * (M * gamma(5 / 2 + n) / R ** 3 / gamma(1 + n)) ** (
            1 / n) / R ** 2  # ~ 2.01
    # For lambda calculation, https://pmocz.github.io/manuscripts/pmocz_sph.pdf

    pos = 2 * np.random.randn(N, 3)  # initial random position for each particle
    vel = np.zeros_like(pos)
    acc = cal_acc(pos, vel, h, Lambda, nu, m, k, n)

    points = [vp.sphere(pos=vp.vec(*point_pos), radius=r) for point_pos in pos]
    color_change = 20
    i = 0
    while True:
        vp.rate(120)

        # Leapfrog integration
        vel += acc * dt / 2
        pos += vel * dt
        acc = cal_acc(pos, vel, h, Lambda, nu, m, k, n)
        vel += acc * dt / 2

        # update particle position for rendering
        for particle_pos, particle in zip(pos, points):
            particle.pos = vp.vec(*particle_pos)

        # update particle's color representing density. High density region-> yellow, low density region -> red
        if i % color_change == 0:
            # Set color according to density
            density = cal_density(pos, pos, m, h)[:, 0]  # [N]
            norm = matplotlib.colors.Normalize(vmin=min(density), vmax=max(density), clip=True)
            mapper = cm.ScalarMappable(norm=norm, cmap=cm.autumn)
            new_colors = [mapper.to_rgba(rho)[:3] for rho in density]
            for particle_color, particle in zip(new_colors, points):
                particle.color = vp.vec(*particle_color)
            i = 0
        i += 1


if __name__ == '__main__':
    main()
