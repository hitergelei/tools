#!/usr/bin/env python
import numpy as np

# BN, AlN, GaN default Params
titles     = ['h-BN'    , 'h-AlN'   , 'h-GaN'   ]
a          = [  2.517   ,   3.113   ,   3.243   ]
m_A        = [ 10.81    ,  26.98    ,  69.72    ]
m_B        = [ 14.01    ,  14.01    ,  14.01    ]
A_0_AA_in  = [ 49.27    ,  32.25    ,  29.56    ]
A_0_AA_out = [ 12.11    ,   5.342   ,   3.663   ]
A_0_BB_in  = [ 53.27    ,  32.14    ,  30.44    ]
A_0_BB_out = [  7.701   ,   2.720   ,   1.393   ]
A_1_AB_LL  = [-21.24    , -16.64    , -15.86    ]
A_1_AB_TT  = [ -7.904   ,  -2.521   ,  -1.728   ]
A_1_AB_ZZ  = [ -4.954   ,  -1.914   ,  -1.051   ]
A_2_AA_LL  = [ -4.664   ,  -2.488   ,  -2.542   ]
A_2_AA_TT  = [  2.814   ,   1.205   ,   1.428   ]
A_2_AA_ZZ  = [  0.4579  ,   0.02102 ,  -0.1795  ]
A_2_AA_LT  = [ -0.9175  ,  -0.2863  ,  -0.5366  ]
A_2_BB_LL  = [ -4.912   ,  -2.045   ,  -2.0155  ]
A_2_BB_TT  = [  1.729   ,   0.7054  ,   0.4856  ]
A_2_BB_ZZ  = [  1.193   ,   0.5038  ,   0.2520  ]
A_2_BB_LT  = [ -0.2871  ,  -0.1949  ,  -0.1133  ]
# A_2_AA_TT  = [0., 0., 0.]
# A_2_AA_LL  = [0., 0., 0.]
# A_2_AA_ZZ  = [0., 0., 0.]
# A_2_AA_LT  = [0., 0., 0.]
# A_2_BB_TT  = [0., 0., 0.]
# A_2_BB_LL  = [0., 0., 0.]
# A_2_BB_ZZ  = [0., 0., 0.]
# A_2_BB_LT  = [0., 0., 0.]
zz_factor = [1.0, 1.0]
# zz_factor = [2.00, 0.5]
mass_factor = None
# mass_factor = 5.0
ylim = ((-7, 47), (-5, 31), (-3, 28))

def get_D_q(
    qa1, qa2,
    m_A, m_B,
    A_0_AA_in, A_0_AA_out,
    A_0_BB_in, A_0_BB_out,
    A_1_AB_TT, A_1_AB_LL, A_1_AB_ZZ,
    A_2_AA_TT, A_2_AA_LL, A_2_AA_ZZ, A_2_AA_LT,
    A_2_BB_TT, A_2_BB_LL, A_2_BB_ZZ, A_2_BB_LT,
    ):
    """
    qa1 = np.dot(q, a1)
    """

    #
    t1 = 2*np.cos(qa1) + np.cos(np.pi/3)**2 *(np.cos(qa1+qa2) + np.cos(qa2)) + 0j
    t2 = 2*np.sin(np.pi/3)**2 *(np.cos(qa1+qa2) + np.cos(qa2)) + 0j
    t3 = 2*np.sin(np.pi/3)*np.cos(np.pi/3)*(np.cos(qa1+qa2) - np.cos(qa2)) + 0j
    t4 = 2j*(np.sin(qa1) - np.sin(qa1+qa2) + np.sin(qa2)) + 0j
    t5 = 2*(np.cos(qa1) + np.cos(qa1+qa2) + np.cos(qa2)) + 0j
    #
    D_q_AAXX = 1/m_A *(A_0_AA_in + t1*A_2_AA_LL + t2*A_2_AA_TT)
    D_q_AAYY = 1/m_A *(A_0_AA_in + t2*A_2_AA_LL + t1*A_2_AA_TT)
    D_q_AAXY = 1/m_A *(t3*(A_2_AA_LL - A_2_AA_TT) + t4*A_2_AA_LT)
    D_q_AAZZ = 1/m_A *(A_0_AA_out + t5*A_2_AA_ZZ)
    #
    D_q_BBXX = 1/m_B *(A_0_BB_in + t1*A_2_BB_LL + t2*A_2_BB_TT)
    D_q_BBYY = 1/m_B *(A_0_BB_in + t2*A_2_BB_LL + t1*A_2_BB_TT)
    D_q_BBXY = 1/m_B *(t3*(A_2_BB_LL - A_2_BB_TT) + t4*A_2_BB_LT)
    D_q_BBZZ = 1/m_B *(A_0_BB_out + t5*A_2_BB_ZZ)
    #
    p1 = np.exp(1j/3*(qa1-qa2)) / np.sqrt(m_A * m_B)
    t6 = (1 + np.exp(-1j*qa1)) *np.cos(np.pi/6)**2
    t7 = (1 + np.exp(-1j*qa1)) *np.sin(np.pi/6)**2 + np.exp(1j*qa2)
    t8 = np.cos(np.pi/6)*np.sin(np.pi/6)*(-1 + np.exp(-1j*qa1))
    t9 = 1 + np.exp(1j*qa2) + np.exp(-1j*qa1)
    #
    D_q_ABXX = p1 * (t6*A_1_AB_LL + t7*A_1_AB_TT)
    D_q_ABYY = p1 * (t7*A_1_AB_LL + t6*A_1_AB_TT)
    D_q_ABXY = p1 * t8*(A_1_AB_LL - A_1_AB_TT)
    D_q_ABZZ = p1 * t9*A_1_AB_ZZ

    #
    D_q = np.zeros((6,6), dtype='complex')

    # 
    D_q[0,0] = D_q_AAXX
    D_q[1,1] = D_q_AAYY
    D_q[0,1] = D_q_AAXY
    D_q[2,2] = D_q_AAZZ

    #
    D_q[3,3] = D_q_BBXX
    D_q[4,4] = D_q_BBYY
    D_q[3,4] = D_q_BBXY
    D_q[5,5] = D_q_BBZZ

    #
    D_q[0,3] = D_q_ABXX
    D_q[1,4] = D_q_ABYY
    D_q[2,5] = D_q_ABZZ
    D_q[0,4] = D_q_ABXY
    D_q[1,3] = D_q_ABXY

    # Hermitian condition
    for i in range(6):
        for j in range(6):
            if i > j:
                D_q[i,j] = np.conj(D_q[j,i])

    return D_q

def apply_rot_sym(D_q):
    R = np.array([[np.cos(2/3*np.pi), -np.sin(2/3*np.pi), 0],
                  [np.sin(2/3*np.pi),  np.cos(2/3*np.pi), 0],
                  [                0,                  0, 1]])
    R66 = np.zeros((6,6), dtype='float')
    R66[0:3, 0:3] = R
    R66[0:3, 3:6] = R
    R66[3:6, 0:3] = R
    R66[3:6, 3:6] = R
    D_q_set = [D_q]
    for i in range(2):
        D_q_set.append(np.matmul(
            np.matmul(R66, D_q_set[-1]),
            R66.T,
            ))
    return np.mean(D_q_set, axis=0)

def main(
    a,
    m_A, m_B,
    A_0_AA_in, A_0_AA_out,
    A_0_BB_in, A_0_BB_out,
    A_1_AB_TT, A_1_AB_LL, A_1_AB_ZZ,
    A_2_AA_TT, A_2_AA_LL, A_2_AA_ZZ, A_2_AA_LT,
    A_2_BB_TT, A_2_BB_LL, A_2_BB_ZZ, A_2_BB_LT,
    zz_factor = [1.0, 1.0],
    mass_factor = None,
    ylim = None,
    title = None,
    ):

    # Engineer
    A_0_AA_out *= zz_factor[0]
    A_0_BB_out *= zz_factor[1]
    A_1_AB_ZZ *= np.sqrt(zz_factor[0] *zz_factor[1])
    A_2_AA_ZZ *= zz_factor[0]
    A_2_BB_ZZ *= zz_factor[1]
    if mass_factor is not None:
        # m_A = m_B * mass_factor
        m_B = m_A * mass_factor

    #
    a1 = np.array([a                  , 0.                 ])
    a2 = np.array([a*np.cos(2/3*np.pi), a*np.sin(2/3*np.pi)])

    R = np.array([[0, -1], [1, 0]])
    b1 = 2*np.pi / np.dot(a1, np.matmul(R, a2)) * np.matmul(R, a2)
    b2 = 2*np.pi / np.dot(a2, np.matmul(R, a1)) * np.matmul(R, a1)

    #
    from ase.dft.kpoints import ibz_points
    points = ibz_points['hexagonal2']
    G = points['Gamma'][:2]
    M = points['M'][:2]
    K = points['K'][:2]
    kpath = [[G, K], [K, M], [M, G]]
    k_labels = ['$\Gamma$', 'K', 'M', '$\Gamma$']
    #
    n_q = 100
    from ss_phonopy import make_band
    # kpoints.shape = (len(kpath), n_q+1, 2)
    kpoints = np.array(make_band(kpath, n_q))

    # Length
    lengths = []
    for i in range(len(kpath)):
        c = np.array(kpath[i][1]) - np.array(kpath[i][0])
        lengths.append(np.linalg.norm(c[0] * b1 + c[1] * b2))

    k_line = []
    for i in range(len(lengths)):
        k_line.append(np.arange(n_q+1) /n_q *lengths[i] + np.sum(lengths[:i]))
    k_line = np.array(k_line)

    #
    w2 = []
    eig_vec = []
    for i in range(len(kpath)):
        _w2 = []
        _ev = []
        for j in range(n_q+1):
            D_q = get_D_q(
                2*np.pi*kpoints[i,j,0],
                2*np.pi*kpoints[i,j,1],
                m_A, m_B,
                A_0_AA_in, A_0_AA_out,
                A_0_BB_in, A_0_BB_out,
                A_1_AB_TT, A_1_AB_LL, A_1_AB_ZZ,
                A_2_AA_TT, A_2_AA_LL, A_2_AA_ZZ, A_2_AA_LT,
                A_2_BB_TT, A_2_BB_LL, A_2_BB_ZZ, A_2_BB_LT,
                )
            # D_q = apply_rot_sym(D_q) # Doesn't work well.
            __w2, __ev = np.linalg.eigh(D_q)
            _w2.append(__w2)
            _ev.append(__ev)
            # if i == 0 and j == 0:
                # np.save('D_q_at_gamma.npy', D_q)
        w2.append(_w2)
        eig_vec.append(_ev)

    # print(D_q)
    # print(D_q[0:3, 0:3] + D_q[0:3, 3:6] + D_q[3:6, 0:3] + D_q[3:6, 3:6]) # All components are better to be small.

    # w.shape = (len(kpath), n_q, 6)
    from phonopy.units import VaspToTHz
    w = np.sqrt(w2, dtype='complex') *VaspToTHz
    w_r = np.real(w)
    w_i = np.imag(w)
    w = w_r - w_i
    # w.shape = (len(kpath), n_q, 6, 6)
    eig_vec = np.array(eig_vec)

    #
    from matplotlib import pyplot as plt
    plt.figure()
    w_concat = w.reshape((-1,6))
    for i in range(6):
        plt.plot(np.concatenate(k_line), w_concat[:,i], c='k')
    if title is not None:
        plt.title(title, fontsize='x-large')
    # plt.xlabel('', fontsize='x-large')
    plt.ylabel('Frequency (THz)', fontsize='x-large')
    # plt.legend(fontsize='large').set_draggable(True)
    if ylim is not None:
        plt.ylim(ylim)
    plt.xlim((0, k_line[-1,-1]))
    plt.subplots_adjust(left=0.28, bottom=0.35, right=0.70, top=0.80, wspace=0.1, hspace=0.2)
    # plt.yticks(range(0,50,10))
    plt.yticks(range(0,50,10))
    plt.xticks([0.]+k_line[:,-1].tolist(), k_labels)
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.axhline(0, ls=':', c='k', lw=0.5)
    plt.grid(alpha=0.2)
    return plt

for i in range(3):
    plt = main(
        a[i],
        m_A[i], m_B[i],
        A_0_AA_in[i], A_0_AA_out[i],
        A_0_BB_in[i], A_0_BB_out[i],
        A_1_AB_TT[i], A_1_AB_LL[i], A_1_AB_ZZ[i],
        A_2_AA_TT[i], A_2_AA_LL[i], A_2_AA_ZZ[i], A_2_AA_LT[i],
        A_2_BB_TT[i], A_2_BB_LL[i], A_2_BB_ZZ[i], A_2_BB_LT[i],
        zz_factor,
        mass_factor,
        ylim[i],
        titles[i],
        )

plt.figure()
# plt.plot(np.array(a         )/a         [0], label='a'         )
# plt.plot(np.array(m_A       )/m_A       [0], label='m_A'       )
# plt.plot(np.array(m_B       )/m_B       [0], label='m_B'       )
plt.plot(np.array(A_0_AA_in )/A_0_AA_in [0], label='A_0_AA_in' )
plt.plot(np.array(A_0_AA_out)/A_0_AA_out[0], label='A_0_AA_out', ls='--')
plt.plot(np.array(A_0_BB_in )/A_0_BB_in [0], label='A_0_BB_in')
plt.plot(np.array(A_0_BB_out)/A_0_BB_out[0], label='A_0_BB_out', ls='--')
plt.plot(np.array(A_1_AB_TT )/A_1_AB_TT [0], label='A_1_AB_TT')
plt.plot(np.array(A_1_AB_LL )/A_1_AB_LL [0], label='A_1_AB_LL')
plt.plot(np.array(A_1_AB_ZZ )/A_1_AB_ZZ [0], label='A_1_AB_ZZ', ls='--')
plt.legend(fontsize='large').set_draggable(True)
# plt.subplots_adjust(left=0.20, bottom=0.35, right=0.70, top=0.90, wspace=0.1, hspace=0.2)
plt.xticks(range(len(titles)), titles)
plt.tick_params(axis="both",direction="in", labelsize='x-large')
plt.grid(alpha=0.5)

plt.figure()
plt.plot(np.array(A_2_AA_TT )/A_2_AA_TT [0], label='A_2_AA_TT')
plt.plot(np.array(A_2_AA_LL )/A_2_AA_LL [0], label='A_2_AA_LL')
plt.plot(np.array(A_2_AA_ZZ )/A_2_AA_ZZ [0], label='A_2_AA_ZZ', ls='--')
plt.plot(np.array(A_2_AA_LT )/A_2_AA_LT [0], label='A_2_AA_LT')
plt.plot(np.array(A_2_BB_TT )/A_2_BB_TT [0], label='A_2_BB_TT')
plt.plot(np.array(A_2_BB_LL )/A_2_BB_LL [0], label='A_2_BB_LL')
plt.plot(np.array(A_2_BB_ZZ )/A_2_BB_ZZ [0], label='A_2_BB_ZZ', ls='--')
plt.plot(np.array(A_2_BB_LT )/A_2_BB_LT [0], label='A_2_BB_LT')
plt.legend(fontsize='large').set_draggable(True)
# plt.subplots_adjust(left=0.28, bottom=0.35, right=0.70, top=0.80, wspace=0.1, hspace=0.2)
plt.xticks(range(len(titles)), titles)
plt.tick_params(axis="both",direction="in", labelsize='x-large')
plt.grid(alpha=0.5)
plt.show()
