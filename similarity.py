#!/usr/bin/env python

## Global modules
import numpy as np
from ase.io import read
from subprocess import call


## Hyper params
load = True
# line_up = False
line_up = True
# Fgpt params
predict_file_list = [
    'GST225-5_invfer-100K.traj',
    # 'GST225-1_kooi-100K-small.traj',
    ]
train_file_list = [
    # 'GST225-1_kooi-100K-small.traj',
    # 'GST225-2_Petrov-100K-small.traj',
    '7.7kooiNferGeTe-w-4phase.traj',
    'rlxed-ran.traj',
    ]
multipole_order = [2]
cutoff_radi     = 6.0
cutoff_num      = [14, 14, 35]

# Histogram
# x_lim = None
x_lim = (0.173412155158, 0.173412155160)
# y_lim_upper = None
y_lim_upper = (0.,400.)
y_lim_lower = (0.,15.)
bins  = 100

# Reliability check params
object_spec   = 'Ge'
if object_spec == 'Ge' or object_spec =='Sb':
    compare_num = 1
elif object_spec == 'Te':
    compare_num = 1
# if object_spec == 'Sb' or object_spec =='Ge':
    # predict_sample_rate = [
        # # 4,
        # 1,
        # ]
    # train_sample_rate   = [
        # # 4,
        # # 4,
        # # 40,
        # 1,
        # 1,
        # 1,
        # ]
# elif object_spec == 'Te':
    # predict_sample_rate = [
        # # 10,
        # 1,
        # ]
    # train_sample_rate   = [
        # # 10,
        # # 10,
        # # 100,
        # 1,
        # 1,
        # 1,
        # ]
# predict_labels = [
    # 'InvFer',
    # # 'Kooi',
    # ]
# train_labels = [
    # 'Kooi',
    # 'Petrov',
    # 'Random',
    # ]

## Preprocess
# File name
p_l_n = ''
for l in predict_file_list:
    p_l_n += l
t_l_n = ''
for l in train_file_list:
    t_l_n += l
m_n = ''
for l in multipole_order:
    m_n += str(l)
c_n = ''
for l in cutoff_num:
    c_n += str(l)
# p_s_r = ''
# for l in predict_sample_rate:
    # p_s_r += str(l)
# t_s_r = ''
# for l in train_sample_rate:
    # t_s_r += str(l)
file_name = 'P-{}_T-{}_mo{}_cr{:.1f}_cn{}_os{}_cn{}_lu{}.npz'.format(p_l_n, t_l_n, m_n, cutoff_radi, c_n, object_spec, compare_num, line_up)#, p_s_r, t_s_r)

## Load
if load == True:
    try:
        # npz = np.load('saved-coverage/'+file_name)
        npz = np.load('saved-similarity/'+file_name)
    except:
        load = False
    else:
        file_name    = npz['file_name']
        # coverage     = npz['coverage']
        sim_avg      = npz['sim_avg']
        line_up      = npz['line_up']

        train_fgpts        = npz['train_fgpts']
        train_Euler_angles = npz['train_Euler_angles']
        # train_labels       = npz['train_labels']
        train_file_list    = npz['train_file_list']
        train_nums         = npz['train_nums']
        train_indices      = npz['train_indices']

        predict_fgpts        = npz['predict_fgpts']
        predict_Euler_angles = npz['predict_Euler_angles']
        # predict_labels       = npz['predict_labels']
        predict_file_list    = npz['predict_file_list']
        predict_nums         = npz['predict_nums']
        predict_indices      = npz['predict_indices']

def get_fgpts(
    cutoff_radi,
    cutoff_num,
    multipole_order,
    file_list,
    # sample_rate,
    ):
    ## Define vector class
    if line_up:
        from dscrpt_lineup import vector
    else:
        from dscrpt_ran import vector
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
        # alist = alist[::sample_rate[i]]
        fgpts_tmp, Euler_tmp = vec.gen_fgpts('alist', alist, rotational_variation=True)
        fgpts.extend(fgpts_tmp.reshape((-1,fgpt_vec_size)))
        Euler_angles.extend(Euler_tmp.reshape((-1,3)))
        nums.append(int(len(alist)*len(alist[0][np.array(alist[0].get_chemical_symbols())==object_spec])))
        file_ind.extend([file_list[i]] * int(len(alist)*len(alist[0])))
        img_ind.extend(list(np.repeat((np.arange(len(alist))).astype('str'),len(alist[0]))))
        spec_ind.extend(list(alist[0].get_chemical_symbols()       )*len(alist))
        atom_ind.extend(list(np.arange(len(alist[0])).astype('str'))*len(alist))
    
    select = np.array(spec_ind) == object_spec
    fgpts = np.array(fgpts)[select]
    Euler_angles = np.array(Euler_angles)[select]
    indices = np.transpose(np.array([file_ind, img_ind, spec_ind, atom_ind]))[select]
    return fgpts, Euler_angles, nums, indices

## Calculate
if load == False:
    ## Get fgpts
    predict_fgpts, predict_Euler_angles, predict_nums, predict_indices = get_fgpts(
        cutoff_radi,
        cutoff_num,
        multipole_order,
        predict_file_list,
        # predict_sample_rate,
        )
    train_fgpts, train_Euler_angles, train_nums, train_indices = get_fgpts(
        cutoff_radi,
        cutoff_num,
        multipole_order,
        train_file_list,
        # train_sample_rate,
        )

    #### Measure
    # coverage = []
    # for i in range(len(predict_fgpts)):
        # dist_norm = np.linalg.norm(predict_fgpts[i] - train_fgpts, axis=-1) / np.linalg.norm(predict_fgpts[i])
        # coverage.append(np.mean(np.sort(dist_norm)[:compare_num]))
    sim_avg = []
    print(predict_fgpts.shape, train_fgpts.shape)
    for i in range(len(predict_fgpts)):
        similarity = np.exp(-np.linalg.norm(predict_fgpts[i] - train_fgpts, axis=-1))
        # Normalize
        # Not implemented yet.
        sim_avg.append(np.mean(np.sort(similarity)[-compare_num:]))

    ## Save
    # call(['mkdir saved-coverage'], shell=True)
    call(['mkdir saved-similarity'], shell=True)
    np.savez(
        # 'saved/'+file_name,
        'saved-similarity/'+file_name,

        file_name = file_name,
        # coverage  = coverage,
        sim_avg   = sim_avg,
        line_up   = line_up,

        train_fgpts        = train_fgpts,
        train_Euler_angles = train_Euler_angles,
        # train_labels       = train_labels,
        train_file_list    = train_file_list,
        train_nums         = train_nums,
        train_indices      = train_indices,

        predict_fgpts        = predict_fgpts,
        predict_Euler_angles = predict_Euler_angles,
        # predict_labels       = predict_labels,
        predict_file_list    = predict_file_list,
        predict_nums         = predict_nums,
        predict_indices      = predict_indices,
        )

## Postprocess
# avrg = np.mean(coverage)
# std  = np.std(coverage)
avrg = np.mean(sim_avg)
std  = np.std(sim_avg)
print('P={} T={} multipole_order={} cutoff_radi={:.1f} cutoff_num={} object_spec={} compare_num={} line_up={}'.format(p_l_n, t_l_n, m_n, cutoff_radi, c_n, object_spec, compare_num, line_up))#, p_s_r, t_s_r))
# print('coverage = {:6.4f} (Lower is similar)'.format(avrg))
print('Average similarity = {:6.4f} (Higher is similar)'.format(avrg))

## Plot
import matplotlib.pyplot as plt
font = {'family':'Arial'}
plt.rc('font', **font)
fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)
# ax: Subplot only for ylabel
ax = fig.add_subplot(111, frameon=False)
ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
ax.grid(False)

## Histogram
# n, bins, patches = ax1.hist(coverage, bins=bins, range=x_lim, facecolor='purple', alpha=0.70)
# n, bins, patches = ax2.hist(coverage, bins=bins, range=x_lim, facecolor='purple', alpha=0.70)
n, bins, patches = ax1.hist(sim_avg, bins=bins, range=x_lim, facecolor='gray')#, alpha=0.70)
n, bins, patches = ax2.hist(sim_avg, bins=bins, range=x_lim, facecolor='gray')#, alpha=0.70)
max_height = np.sort(n)[-1]
# Define title
title = ''
for i in range(len(predict_file_list)):
    title += 'P--> ' + predict_file_list[i] + '\n'
for i in range(len(train_file_list)):
    title += 'T--> ' + train_file_list[i] + '\n'

ax1.set_title(title+'Normalized Distance Minimum Histo.')
# ax2.set_xlabel('Total %d fgpts, average = %.4f, sigma = %.3f\nCompare number = %d, Object species = %s' % (len(coverage), avrg, std, compare_num, object_spec))
ax2.set_xlabel('Total %d fgpts, average = %.4f, sigma = %.3f\nCompare number = %d, Object species = %s' % (len(sim_avg), avrg, std, compare_num, object_spec))
ax.set_ylabel('population', fontsize='large')
# Deviation bar (start from average)
ax1.barh(7., std, height=.8, left=avrg, color='black')
ax2.barh(7., std, height=.8, left=avrg, color='black')
ax1.set_xlim(x_lim)
# y_lim_upper = (100.,max_height + 50)
ax1.set_ylim(y_lim_upper)
ax2.set_ylim(y_lim_lower)

ax1.grid(alpha=0.2)
ax2.grid(alpha=0.2)
ax1.tick_params(axis="both",direction="in", labelsize='large')
ax2.tick_params(axis="both",direction="in", labelsize='large')
plt.subplots_adjust(top=0.6, hspace=0.1, bottom=0.4)
plt.show()
