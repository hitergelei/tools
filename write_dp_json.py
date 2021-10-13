#!/usr/bin/env python
import json
from sys import argv

# argv sequence: batch_size, stop_batch, disp_freq, save_freq, json_name

#
with open('input.json') as rf:
    data = json.load(rf)

#
# data['training']['batch_size'] = int(argv[1])
# data['training']['stop_batch'] = int(argv[2])
data['training']['training_data']['batch_size'] = int(argv[1])
data['training']['numb_steps'] = int(argv[2])
data['training']['disp_freq'] = int(argv[3])
data['training']['save_freq'] = int(argv[4])

#
out_data = json.dumps(data, indent=4)
with open(argv[5], 'w') as wf:
    wf.write(out_data)
