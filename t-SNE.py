#!/usr/bin/env python

## Global modules
import numpy as np
from ase.io import read
from subprocess import call

## Global params
load = True
ticks = False
# Fgpt params
file_list = [
    'GST225-1_kooi-100K-small.traj',
    'GST225-2_Petrov-100K-small.traj',
    'GST225-5_invfer-100K.traj',
    'rlxed-ran.traj',
    ]
labels          = [
    'Kooi',
    'Petrov',
    'InvFer',
    'Random',
    ]
multipole_order = [2]
cutoff_radi     = 6.0
cutoff_num      = [14, 14, 35]
# t-SNE params
object_spec   = 'Ge'
dimension     = 2
if object_spec == 'Sb' or object_spec =='Ge':
    sample_rate   = [
        4,
        4,
        4,
        40,
        ]
elif object_spec == 'Te':
    sample_rate   = [
        10,
        10,
        10,
        100,
        ]
n_iter        = 5000
learning_rate = 50
# perplexity    = 2
perplexity    = 5
# perplexity    = 30
# perplexity    = 50
# perplexity    = 100
# perplexity    = 200

## Preprocess
# File name
l_n = ''
for l in labels:
    l_n += l
m_n = ''
for l in multipole_order:
    m_n += str(l)
c_n = ''
for l in cutoff_num:
    c_n += str(l)
s_n = ''
for l in sample_rate:
    s_n += str(l)
file_name = '{}_mo{}_cr{:.1f}_cn{}_os{}_d{:d}_sr{}_ni{:d}_lr{:d}_p{:d}.npz'.format(l_n, m_n, cutoff_radi, c_n, object_spec, dimension, s_n, n_iter, learning_rate, perplexity)

## Load
if load == True:
    try:
        npz = np.load('saved/'+file_name)
    except:
        load = False
    else:
        file_name    = npz['file_name']
        projection   = npz['projection']
        fgpts        = npz['fgpts']
        Euler_angles = npz['Euler_angles']
        labels       = npz['labels']
        nums         = npz['nums']
        indices      = npz['indices']

## Calculate
if load == False:
    ## Define vector class
    from dscrpt import vector
    vec = vector(
        cutoff_radi     = cutoff_radi,
        cutoff_num      = cutoff_num,
        multipole_order = multipole_order+['r'],
        # one_box       = True,
        # logfile_name  = 'log.txt',
        )

    ## Get fgpts
    nums           = []
    file_ind       = []
    img_ind        = []
    spec_ind       = []
    atom_ind       = []
    fgpts          = []
    Euler_angles   = []
    fgpt_vec_size = np.sum(cutoff_num) * (1+ len(multipole_order)*3)
    for i in range(len(file_list)):
        print("   >> Generating {:}'s fgpts <<".format(file_list[i]))
        alist = read(file_list[i], ':')
        alist = alist[::sample_rate[i]]
        fgpts_tmp, Euler_tmp = vec.gen_fgpts('alist', alist, rotational_variation=True)
        fgpts.extend(fgpts_tmp.reshape((-1,fgpt_vec_size)))
        Euler_angles.extend(Euler_tmp.reshape((-1,3)))
        nums.append(int(len(alist)*len(alist[0][np.array(alist[0].get_chemical_symbols())==object_spec])))
        file_ind.extend([labels[i]] * int(len(alist)*len(alist[0])))
        img_ind.extend(list(np.repeat((np.arange(len(alist))*sample_rate[i]).astype('str'),len(alist[0]))))
        spec_ind.extend(list(alist[0].get_chemical_symbols()       )*len(alist))
        atom_ind.extend(list(np.arange(len(alist[0])).astype('str'))*len(alist))
    
    select = np.array(spec_ind) == object_spec
    fgpts = np.array(fgpts)[select]
    Euler_angles = np.array(Euler_angles)[select]
    indices = np.transpose(np.array([file_ind, img_ind, spec_ind, atom_ind]))[select]

    #### t-SNE
    from sklearn.manifold import TSNE
    model = TSNE(
        n_components  = dimension,
        n_iter        = n_iter,
        learning_rate = learning_rate,
        perplexity    = perplexity,
        verbose       = 1,
        init          = 'pca',
        random_state  = 25563,
        )
    projection = model.fit_transform(fgpts)

    ## Save
    call(['mkdir saved'], shell=True)
    np.savez(
        'saved/'+file_name,
        file_name    = file_name,
        projection   = projection,
        fgpts        = fgpts,
        Euler_angles = Euler_angles,
        labels       = labels,
        nums         = nums,
        indices      = indices,
        )

## Postprocess
x1 = projection[:,0]
x2 = projection[:,1]
if dimension == 3:
    x3 = projection[:,2]

## Plot
import matplotlib.pyplot as plt
font = {'family':'Arial'}
plt.rc('font', **font)
plt.tick_params(axis="y",direction="in")
plt.tick_params(axis="x",direction="in")

count = 0
if dimension == 3:
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(nums)):
        ax.scatter(x1[count:count+nums[i]],x2[count:count+nums[i]],x3[count:count+nums[i]], label=labels[i])
        count += nums[i]
    # ax.legend(scatterpoints = 1, fontsize='large', loc='upper right')
elif dimension == 2:
    for i in range(len(nums)):
        plt.scatter(x1[count:count+nums[i]],x2[count:count+nums[i]], s=20, alpha=1.0, label=labels[i])
        count += nums[i]
    plt.legend(scatterpoints = 1, fontsize='large', loc='upper right')
if not ticks:
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False)
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        right=False,      # ticks along the bottom edge are off
        left=False,         # ticks along the top edge are off
        labelleft=False)

plt.subplots_adjust(left=0.15, bottom=0.05, right=0.85, top=0.90, wspace=0.2, hspace=0.2)
plt.title(object_spec+', Perplexity='+str(perplexity), fontdict={'fontsize':20})
plt.show()
