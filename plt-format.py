
from matplotlib import pyplot as plt
font = {'family':'sans-serif', 'sans-serif':'Arial'}
plt.rc('font', **font)

plt.title('', fontsize='x-large', pad=None)
plt.xlabel('', fontsize='x-large')
plt.ylabel('', fontsize='x-large')
# plt.xscale('log')
plt.tick_params(axis="both",direction="in", labelsize='x-large')
plt.subplots_adjust(left=0.30, bottom=0.30, right=0.70, top=0.70, wspace=0.20, hspace=0.20)
plt.legend(fontsize='large').set_draggable(True)
plt.grid(alpha=0.5)
