import numpy as np

def get_rotational_matrix(
    new_basis,
    ref_basis = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],
    ):
    """

    Formula >>> a' = R @ a
    new_basis (=a') - Target basis set
    ref_basis (=a)  - Initial basis set
    R               - Rotational matrix of coordinate rotation
                       (If ref_basis is an identity, R = a'. (Rotational matrix itself)
    !!! NOTE !!! As you can see this function is almost useless. Just for notation purpose.

e.g.   Ref.             New basis
    [[1,0,0],          [[x1,x2,x3],
     [0,1,0],    to     [y1,y2,y3],
     [0,0,1]]           [z1,z2,z3]]

    """
    R = np.matmul(new_basis, np.linalg.inv(ref_basis))
    return R

def get_Euler_angles(vector):
    """
    vector - 3x3 matrix like described below. (Old basis representation)
    Return "proper" Euler angles (ZXZ rotation) in radian unit.
       Old              New (input)
    [[1,0,0],          [[x1,x2,x3],
     [0,1,0],    to     [y1,y2,y3],
     [0,0,1]]           [z1,z2,z3]]
    https://en.wikipedia.org/wiki/Euler_angles
    Caution) x, y, z each must be unit vector.
    Note!!!! Euler cosines can be multimapped to vectors
    Recommend to use rotation matrix instead
    """
    # Check input
    V = np.array(vector, dtype='float')
    # if np.shape(np.squeeze(V)) != (3,3):
        # raise ValueError('Vector is not in shape of (3,3) but in {}'.format(str(np.shape(V))))
    # Get angles
    a = np.arccos(-V[2,1] / np.sqrt(1. - V[2,2]**2))
    b = np.arccos( V[2,2])
    c = np.arccos( V[1,2] / np.sqrt(1. - V[2,2]**2))
    return np.array([a, b, c])

def Euler_rotation(vectors, angles, inverse=False):
    """
    vectors - List of vectors of shape (3,)
    Return vectors in new bases after 'axis' Euler rotation with given angles.
    Angle must be given in proper Euler angles. (ZXZ)
    Vectors can be given in the form of list or tuple format of shape (n,3) or also single vector of shape (3,).
    Cf) some_ang + inverse=True == -some_ang[::-1] + inverse=False

    >>Caution!!<<
    some_ang  --> vector rotation
    -some_ang --> coordinate rotation
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

