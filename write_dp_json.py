#!/usr/bin/env python
import json
from sys import argv

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
data['training']['validation_data']['batch_size'] = int(argv[5])
data['training']['validation_data']['numb_btch'] = int(argv[6])

#
out_data = json.dumps(data, indent=4)
with open(argv[7], 'w') as wf:
    wf.write(out_data)
