#!/usr/bin/env python

def get_gpu_avail(
    ):
    from subprocess import check_output

    # Get output
    out = str(check_output('nvidia-smi')).split('\\n')

    # Find first line
    line_num = []
    for i in range(len(out)):
        if out[i][:2] == '|=':
            line_num.append(i+1)

    # Get indices of GPUs.
    gpu_total = set()
    while True:
        words = out[line_num[0]].split()
        if len(words) == 0:
            break
        gpu_total.add(words[1])
        line_num[0] += 3

    #Get indices of GPUs in use.
    gpu_in_use = set()
    while True:
        words = out[line_num[1]].split()
        if len(words) == 1:
            break
        if words[1] == 'No':
            break
        gpu_in_use.add(words[1])
        line_num[1] += 1

    return list(gpu_total - gpu_in_use)

# @ Main
gpu_avail = get_gpu_avail()
if len(gpu_avail) == 0:
    raise RuntimeError('No GPU device available.')
print(get_gpu_avail()[0])
