#!/usr/bin/env python
import numpy as np
from subprocess import call

call('rm -rf train test', shell=True)
call('mkdir train test', shell=True)
call('mv set.000 train/', shell=True)
call('mv set.001 test/set.000', shell=True)
call('cp type.raw train/', shell=True)
call('cp type.raw test/', shell=True)
