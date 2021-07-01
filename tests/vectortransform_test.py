from numpy import *
import math
import matplotlib.pyplot as plt


def normvec(vec):
    a = empty(shape(vec))
    for i in range(len(vec)):
        a[i] = unit(vec[i])
    return a


def scal(v1, v2):
    return vdot(v1, v2)


def unit(v):
    return 1 / sqrt(scal(v, v)) * v


def gram_schmidt(L1):
    L2 = [unit(L1[0])]
    for i in range(1, len(L1)):
        s = L1[i]
        for j in range(i):
            s = s - scal(L1[i], L2[j]) * L2[j]
        L2.append(unit(s))
    return array(L2)


def transform_matrix(f1, f2):  # transforms from f1 to f2
    R = array([
        [dot(f2[0], f1[0]), dot(f2[0], f1[1]), dot(f2[0], f1[2])],
        [dot(f2[1], f1[0]), dot(f2[1], f1[1]), dot(f2[1], f1[2])],
        [dot(f2[2], f1[0]), dot(f2[2], f1[1]), dot(f2[2], f1[2])]
    ])
    return R


def rot_vector_3D(_, system, theta, phi):
    v_new_in_system = [math.sin(phi) * math.cos(theta),
                       math.sin(phi) * math.sin(theta),
                       math.cos(phi)]
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    R2 = transform_matrix(system, ECEF)

    rotated_vector = R2.dot(v_new_in_system)
    # print(rotated_vector)
    return rotated_vector


def angle_between(v1, v2):
    v1 = unit(array(v1))
    v2 = unit(array(v2))
    return degrees(round(arccos(round(dot(v1, v2), 5)), 5))


def get_rotated_vectors_from_V_to_earths_system(s1, s2, s3, rot_theta, rot_phi):
    system = array([s1, s2, s3])
    system = flip(gram_schmidt(flip(array(system))))
    # system = gram_schmidt(array(system))

    system = normvec(system)
    v_rotated = rot_vector_3D(None, system, rot_theta, rot_phi)  # v_rotated mar a fold rendszereben van
    # for i in range(len(system)):
    # 	print('base vector {}: '.format(i), angle_between(allo_cs[i],system[i]))

    return v_rotated


def plot_cmap(color_matrix, phi_max, theta_max):
    fig = plt.imshow(color_matrix, extent=[0, degrees(phi_max), 0, degrees(theta_max)], aspect='auto')

    plt.title('24h_average, SD rotated')
    # plt.legend()
    # plt.tight_layout()
    plt.ylabel('Theta')
    plt.xlabel('Phi')
    plt.colorbar(fig)

    plt.show()


allocsillag_r = normvec(random.randint(0, 100, 9).reshape((3, 3)))
# print('a :', allocsillag_r)

sd = allocsillag_r[2]

theta = radians(10)
phi = radians(45)

# sd_rot = get_rotated_vectors_from_V_to_earths_system(sd, allocsillag_r[0], allocsillag_r[1], theta, phi)

# print(angle_between(sd, sd_rot))

# a = array([[1,2,3],[4,5,6],[7,8,9]])
# af = flip(a)

# for i in range(len(a)):
# 	print(angle_between(a[i],af[i]))

# print(a)
# print(flip(a))
# print(flip(flip(a)))

# # print(gram_schmidt(flip(array(a))))

# # print(flip(gram_schmidt(flip(array(a)))))
# =========================================================================================================================================================

# V az allocsillagok rendszere a mi rendszerunkbol
V = gram_schmidt(normvec(array([[1, 0, 2], [0, 3, 5], [9, 4, 0]])))
print('Orto. frame:', V)

print(linalg.norm(V[0]))
print(linalg.norm(V[1]))
print(linalg.norm(V[2]))

# I a keresett irany mindket rendszerben
# I = unit(array([1,1,2]))
I = unit(array([1, 1, 1]))

resolution = radians(5.0)
theta_max = math.pi * 2.0  # + resolution * 0.1
phi_max = math.pi
rot_theta = arange(0.0, theta_max, resolution)
rot_phi = arange(0.0, phi_max, resolution)

cmap = empty((len(rot_theta), len(rot_phi)))
for j in range(len(rot_phi)):
    for i in range(len(rot_theta)):
        theta = rot_theta[i]
        phi = rot_phi[j]
        cmap[i][j] = dot(I, get_rotated_vectors_from_V_to_earths_system(V[0], V[1], V[2], theta, phi))
    # cmap[i][j] = dot(I, get_rotated_vectors_from_V_to_earths_system(V[0], V[1], V[2], theta, phi))
    # cmap[i][j] = scal(V[0], V[2])

plot_cmap(cmap, phi_max, theta_max)
