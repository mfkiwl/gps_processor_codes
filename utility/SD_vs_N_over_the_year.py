from numpy import *
import math



def normvec(vec):
	a = empty(shape(vec))
	for i in range(len(vec)):
		a[i] = unit(vec[i])
	return a

def scal(v1, v2):
	return vdot(v1, v2)

def unit(v):
	return 1/sqrt(scal(v,v))*v


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


def rot_vector_3D(vector, system, theta, phi):
	v_new_in_system = [math.sin(phi) * math.cos(theta),
		        math.sin(phi) * math.sin(theta),
		        math.cos(phi)]
	ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
	R2 = transform_matrix(system, ECEF)
	rotated_vector = R2.dot(v_new_in_system)
	return rotated_vector



def rot_SD_3D(SD_init, s2, s3, rot_theta, rot_phi):
	SD_rotated = empty(shape(SD_init))
	for i in range(len(SD_init)):
		system = flip(gram_schmidt( flip(array([s2[i], s3[i], SD_init[i]]))))
		system = normvec(system)
		SD_rotated[i] = rot_vector_3D(SD_init[i], system, rot_theta, rot_phi)
	return SD_rotated


allocsillag_r = normvec(random.randint(0, 100, 9).reshape((3,3)))
print(allocsillag_r)

sd = allocsillag_r[2]

theta = radians(0)
phi = radians(45)


sd_rot = rot_SD_3D(sd, allocsillag_r[0], allocsillag_r[1], theta, phi)

print(degrees(arccos(dot(sd, sd_rot))))








