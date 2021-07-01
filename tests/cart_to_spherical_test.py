from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import math
from sklearn import preprocessing
import sys

from data_postprocessing.triangle_method_rework.rework_lib import rotate2darray, verify_angle

sys.path.insert(1, '/Users/kelemensz/Documents/Research/GPS/gps_processor_codes/data_postprocessing')


# from triangle_method_rework.rework_lib import get_common, normvec, cart2sph


def scal(v1, v2):
    # return vdot(array(v1), array(v2))
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]


def unit(v):
    return 1.0 / sqrt(scal(v, v)) * array(v)


def normvec(vec):
    a = empty(shape(vec))
    for i in range(len(vec)):
        a[i] = unit(vec[i])
    return a


def cart2sph(x, y, z):
    hxy = hypot(x, y)
    # r = hypot(hxy, z)
    el = arctan2(z, hxy)
    az = arctan2(y, x)
    return az, el


def get_ij_on_map(direction, l_theta, l_phi, resolution):
    theta, phi = cart2sph(direction[0], direction[1], direction[2])
    # print(theta, phi)
    I_f = 0
    J_f = 0
    for i in range(-l_theta, l_theta):
        if i * resolution > theta:
            I_f = i + l_theta
            break
    for i in range(-l_phi, l_phi):
        if i * resolution > phi:
            J_f = i + l_phi
            break
    return I_f, J_f


def build_spherical_cmap_from_3Dvectors(D3vectors, resolution):
    theta_max = math.pi
    phi_max = math.pi / 2.0
    rot_theta = arange(-theta_max, theta_max, resolution)
    rot_phi = arange(-phi_max, phi_max, resolution)
    l_theta = int(len(rot_theta) / 2.0)
    l_phi = int(len(rot_phi) / 2.0)
    cmap_v = zeros((len(rot_theta), len(rot_phi)))
    value = 1
    thetas = []
    phis = []
    for c in D3vectors:
        theta, phi = cart2sph(c[0], c[1], c[2])
        thetas.append(theta)
        phis.append(phi)

        i, j = get_ij_on_map(c, l_theta, l_phi, resolution)
        # print(i,j)
        cmap_v[i][j] += value
    plt.hist(thetas)
    plt.show()
    plt.clf()
    plt.hist(phis, density=True, bins=50)
    # z = np.arange(-pi/2., pi/2., 0.1)
    # plt.plot(z, np.cos(z)/2.0)
    plt.show()
    plt.clf()
    return cmap_v


def plot_cmap(color_matrix, phi_max, theta_max):
    fig = plt.imshow(color_matrix, extent=[0, degrees(theta_max), 0, degrees(phi_max)], aspect='auto')
    # fig = plt.imshow(color_matrix)
    plt.ylabel('Theta')
    plt.xlabel('Phi')
    plt.colorbar(fig)

    plt.show()


def generate_direction_uniformly():
    phi = np.random.uniform(low=-pi, high=pi)
    costheta = np.random.uniform(low=-1.0, high=1.0)
    theta = arccos(costheta)
    # theta = arcsin(costheta)
    r = 1
    return array([r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)])


def generate_directions_uniformly(N):
    # randomvectors = normvec(2*np.random.rand(100000, 3)-1.0)
    # randomvectors = normvec(np.random.randint(-100, 100, size=(10000, 3)))
    # randomvectors = np.random.randint(-100, 100, size=(10000, 3))
    # randomvectors = np.array([[1,0,0],[0,1,0],[0,0,1]])
    randomvectors = np.zeros(shape=[N, 3])
    for i in range(N):
        randomvectors[i] = generate_direction_uniformly()
    return randomvectors


randomvectors = generate_directions_uniformly(100000)

# randomvectors = preprocessing.normalize(randomvectors)

# plt.hist(randomvectors[:, 0])
# plt.show()
# plt.clf()
# plt.hist(randomvectors[:, 1])
# plt.show()
# plt.clf()
# plt.hist(randomvectors[:, 2])
# plt.show()
# plt.clf()

#
cmap = build_spherical_cmap_from_3Dvectors(randomvectors, radians(5))
print(cmap)
#
#
theta_max = math.pi * 2.0  # + resolution * 0.1
phi_max = math.pi

cmap[cmap < 1] = nan
plot_cmap(cmap.T, phi_max, theta_max)

# normed = normvec(randomvectors)
# # verify_angle(randomvectors, normed)
#
# if len(normed) == len(randomvectors):
#     for i in range(len(randomvectors)):
#         # print(angle_between(normed[i], randomvectors[i]))
#         print(np.cross(normed[i], randomvectors[i]))


