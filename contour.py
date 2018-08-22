from pylab import *

nx = 300
ny = 300
data = np.loadtxt('./BERRYCURV.dat')
colorlevel = 51
x,y,z = data[:,0], data[:,1], data[:,3]

xi = np.linspace(x.min(),x.max(),nx)
yi = np.linspace(y.min(),y.max(),ny)
zi = griddata(x,y,z,xi,yi)
maxz=z.max()
minz=z.min()
deltaz=(maxz-minz)/colorlevel
level=arange(minz,maxz+deltaz/5,deltaz)
CS= plt.contourf(xi,yi,zi,level,cmap=cm.RdYlBu)
bounds = np.linspace(z.min(),-1*z.min(),10)
plt.colorbar(ticks=bounds,format='%5.4f')

xin=x.min()
xen=x.max()
yin=y.min()
yen=y.max()
plt.axis([xin,xen,yin,yen])
plt.setp(gca(),yticks=(arange(0)),xticks=(arange(0)),aspect='equal')
setp(gca(), yticklabels=[], yticks=(), xticks=())

savefig('./BERRYCURV.eps',bbox_inches='tight',transparent=True,pad_inches=0)

