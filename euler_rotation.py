import numpy as np

def get_Euler_angles(vectors):
    """
    Return proper Euler angles (ZXZ rotation) in radian unit.
    [[1,0,0],          [[x1,x2,x3],
     [0,1,0],    to     [y1,y2,y3],
     [0,0,1]]           [z1,z2,z3]]
    https://en.wikipedia.org/wiki/Euler_angles
    Caution) x, y, z each must be unit vector.
    """
    # Check input
    V = np.array(vectors, dtype='float')
    # if np.shape(np.squeeze(V)) != (3,3):
        # raise ValueError('Vectors are not in shape of (3,3) but in {}'.format(str(np.shape(V))))
    # Get angles
    a = np.arccos(-V[2,1] / np.sqrt(1. - V[2,2]**2))
    b = np.arccos(V[2,2])
    c = np.arccos(V[1,2] / np.sqrt(1. - V[2,2]**2))
    return np.array([a, b, c])

def Euler_rotation(vectors, angles, inverse=False):
    """
    Return vectors in new bases after 'axis' Euler rotation with given angles.
    Angle must be given in proper Euler angles. (ZXZ)
    Vectors can be given in the form of list or tuple format of shape (n,3) or also single vector of shape (3,).
    """
    # Check input
    V_list = np.array(vectors, dtype='float')
    # if np.shape(np.squeeze(V_list[0])) != (3,):
        # raise ValueError('Vector is not in shape of (3,) but in {}'.format(str(np.shape(V))))
    A = np.array([angles], dtype='float')
    # if np.shape(np.squeeze(A)) != (3,):
        # raise ValueError('Angles are not in shape of (3,) but in {}'.format(str(np.shape(A))))
    # Define Euler rotation
    from scipy.spatial.transform import Rotation as R
    rot = R.from_euler('zxz', A)
    # Get vectors
    return(np.array(rot.apply(V_list, inverse)))

