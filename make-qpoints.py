import numpy as np

# Params
n_mesh = 10
origin = np.array([0., 0., 0.])
b1 = np.array([ 0.37, -0.37,     0])
b2 = np.array([ 0   ,  0.37, -0.37])
# origin = np.array([0.5, 0.5, 0.5])
# b1 = np.array([ 0.75995, 0.24005, 0.50000]) - origin
# b2 = np.array([ 0.50000, 0.75995, 0.24005]) - origin

# First quadrant
q1 = []
for i in range(1,n_mesh+1):
    for j in range(0,n_mesh+1):
        q1.append(i*b1 + j*b2)
q1 = np.array(q1)

# Second quadrant
q2 = []
for i in range(0,n_mesh+1):
    for j in range(0,n_mesh+1):
        if i+j <= n_mesh:
            q2.append(-i*b1 + j*b2)
q2 = np.array(q2)

q = np.concatenate([[[0.,0.,0.]], q1, q2, -q1, -q2], axis=0) /n_mesh
q += origin

np.save('qpoints1.npy', q)
# np.save('qpoints2.npy', q)
