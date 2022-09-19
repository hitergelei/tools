#!/usr/bin/env python
import numpy as np
from subprocess import call

class Voronoi():
    def __init__(
        self,
        atoms_file,
        ):
        """
        atoms_file (ASE readable structure)
            - Only a single structure is read. If file include multiple images, it will load the final structure only.
        """
        import datetime
        now = datetime.datetime.now()
        time = now.strftime('%Y-%m-%d %H:%M:%S')
        print("""
        ### Reference: 
        ### Voronoi-volume-matrix code by Young-Jae Choi @ POSTECH, Republic of Korea
        ### Code start time: {}
        """.format(time))

        self.atoms_file = atoms_file
        from ase.io import read
        self.atoms = read(atoms_file, -1)
        self.dynmat = None
        self.w2 = None
        self.eps = None
        self.f_box = None
        self.f_voro = None
        self.neigh_list = None
        self.facet_list = None
        self.A = None
        self.ul = None

    def set_dynmat(self, dynmat):       
        """
        """
        self.dynmat = dynmat

    def read_fix_phonon_log(
        self,
        log_file,
        ):
        with open(log_file) as f:
            lines = f.readlines()
        raw = np.array(lines[-1][:-1].split(), dtype=float)[3:] #*1e-6 /8.61734e-05
        dynmat = raw[::2] + 1.j*raw[1::2]
        N = int(np.sqrt(dynmat.shape[0])/3.)
        # dynmat = np.transpose(dynmat.reshape(N, N, 3, 3), (0,2,1,3)).reshape(3*N, 3*N)
        dynmat = dynmat.reshape(3*N, 3*N)
        self.dynmat = dynmat

    def set_eigen_sets(self, w2, eps):
        self.w2 = w2
        self.eps = eps

    def get_eigen_sets(self):
        if self.w2 == None and self.eps == None:
            if self.dynmat == None:
                raise RuntimeError('Either eigen set or dynamical matrix has to be set in advance.')
            else:
                self.w2, self.eps = np.linalg.eig(self.dynmat)
        return self.w2, self.eps

    def compute_voronoi(self):
        if self.f_voro == None:
            from freud.box import Box
            f_box = Box.from_matrix(np.array(self.atoms.get_cell()))
            self.f_box = f_box
            from freud.locality import Voronoi as freud_Voronoi
            f_voro = freud_Voronoi()
            f_voro.compute((f_box, self.atoms.get_positions()))
            self.f_voro = f_voro
            self.polytopes_wrap = []
            for i in range(len(f_voro.polytopes)):
                self.polytopes_wrap.append(self.f_box.wrap(self.f_voro.polytopes[i]))

    def get_neigh_list(self):
        if self.neigh_list is None:
            self.compute_voronoi()
            neigh_list = []
            neigh_dist = []
            s_ij = []
            for i in range(len(self.f_voro.nlist.segments)):
                neigh_list.append(self.f_voro.nlist[self.f_voro.nlist.segments[i]:self.f_voro.nlist.segments[i]+self.f_voro.nlist.neighbor_counts[i], 1].tolist())
                neigh_dist.append(self.f_voro.nlist.distances[self.f_voro.nlist.segments[i]:self.f_voro.nlist.segments[i]+self.f_voro.nlist.neighbor_counts[i]])
                s_ij.append(self.f_voro.nlist.weights[self.f_voro.nlist.segments[i]:self.f_voro.nlist.segments[i]+self.f_voro.nlist.neighbor_counts[i]])
            self.neigh_list = neigh_list
            self.neigh_dist = neigh_dist
            self.s_ij = s_ij

    def get_facet_list(self):
        if self.facet_list is None:
            self.get_neigh_list()
            print('Calculating facet_list', flush=True)
            from shapely.geometry import Polygon
            from scipy.spatial import ConvexHull
            facet_list = []
            centroid_list = []
            for i in range(len(self.neigh_list)):
                print(' ({:5.1f}%)  facet_list ind: {}'.format((i+1)/len(self.neigh_list)*100., i+1), flush=True)
                facet_list.append([])
                centroid_list.append(np.zeros((len(self.neigh_list[i]),3)))
                for j in range(len(self.neigh_list[i])):
                    facet_list[i].append([])
                    for k in range(len(self.f_voro.polytopes[i])):
                        if np.amin(np.linalg.norm(self.polytopes_wrap[self.neigh_list[i][j]] - self.polytopes_wrap[i][k], axis=-1)) < 1e-4:
                            facet_list[i][j].append(self.f_voro.polytopes[i][k])
                    facet_list[i][j] = np.array(facet_list[i][j])
                    try:
                        facet_list[i][j] = facet_list[i][j][ConvexHull(facet_list[i][j][:,0:2]).vertices]
                    except:
                        facet_list[i][j] = facet_list[i][j][ConvexHull(facet_list[i][j][:,1:3]).vertices]
                    centroid_list[i][j][0:2] = list(Polygon(facet_list[i][j][:, 0:2]).centroid.coords)[0]
                    centroid_list[i][j][1:3] = list(Polygon(facet_list[i][j][:, 1:3]).centroid.coords)[0]
            self.facet_list = facet_list
            self.centroid_list = centroid_list
            print('facet_list calculation completed', flush=True)

    def get_A_mat(self, npy_name='A_mat.npy'):
        r"""

                                    S_{ij}
        A_{i,j\alpha} (i != j) = ------------ * (r_{j\alpha} - c_{ij\alpha})
                                  V_i r_{ij}

        A_{i,j\alpha} (i == j) = - (sum)_{j != i} A_{i,j\alpha}

        where c_{ij\alpha} is the centroid of the facet polygon S_{ij}.

        """
        if self.A is None:
            try:
                A = np.load(npy_name)
            except:
                do_calc = True
            else:
                do_calc = False

            if do_calc:
                self.get_facet_list()
                posi = self.atoms.get_positions()
                A = np.zeros((len(self.atoms), len(self.atoms), 3))
                for i in range(len(self.neigh_list)):
                    for j in range(len(self.neigh_list[i])):
                        nei_i = self.neigh_list[i][j]
                        A[i,nei_i] = self.s_ij[i][j] /self.f_voro.volumes[i] /self.neigh_dist[i][j] *self.f_box.wrap(posi[nei_i]-self.centroid_list[i][j])
                    A[i,i] = -np.sum(A[i], axis=0)
                A = A.reshape((len(self.atoms), 3*len(self.atoms)))
                np.save(npy_name, A)
            self.A = A

    def LT_decompose(
        self,
        ):
        """
        u_L = A^T (A A^T)^{-1} A u
        u_T = 1 - u_L
        """

        if self.ul is None:
            A = self.A
            ul = A.T @ np.linalg.inv(A @A.T) @A @self.eps.T
            ut = self.eps.T -ul
            self.ul = ul.T
            self.ut = ut.T

    def plot_LT_VDOS(
        self,
        nbins=200,
        nan_to_zero=True,
        plot_imaginary=True,
        gsmear_std=0.,
        xlim=(1e-5, None),
        xlim_bp=(1e-5, None),
        ylim=(None, None),
        ):
        """
        gsmear_std (float)
            - Sigma (std) for Gaussian smearing in THz unit.
        """

        self.LT_decompose()
        w = np.sqrt(self.w2)
        from ase.units import fs
        w *= (fs*1e+3/2/np.pi)

        if nan_to_zero:
            w = np.where(np.isnan(w), 0, w)
        eigval_plot = np.real(w)
        if plot_imaginary:
            eigval_plot -= np.where(np.imag(w) < 0, np.imag(w), 0)

        hist_l, edges = np.histogram(eigval_plot, bins=nbins, weights=np.sum(self.ul**2, axis=-1))
        hist_t, edges = np.histogram(eigval_plot, bins=nbins, weights=np.sum(self.ut**2, axis=-1))
        hist_tot, edges = np.histogram(eigval_plot, bins=nbins, density=True)
        mids = (edges[0:-1] + edges[1:]) /2.
        dw = mids[1] -mids[0]
        hist_l /= len(w) *dw
        hist_t /= len(w) *dw

        if gsmear_std not in (0., False, None):
            print(' Gaussian smearing...', flush=True)
            from scipy.ndimage.filters import gaussian_filter1d
            hist_l = gaussian_filter1d(hist_l, gsmear_std /dw)
            hist_t = gaussian_filter1d(hist_t, gsmear_std /dw)
            hist_tot = gaussian_filter1d(hist_tot, gsmear_std /dw)

        #
        from matplotlib import pyplot as plt
        for i in range(2):
            fig, ax = plt.subplots()
            if i==0:
                plt.plot(hist_l, mids, c='r', label='L')
                plt.plot(hist_t, mids, c='b', label='T')
                plt.plot(hist_tot, mids, c='k')
                plt.xlabel('$g(\omega)$', fontsize='x-large')
                plt.xlim(xlim)
            else:
                plt.plot(hist_l/mids**2, mids, c='r', label='L')
                plt.plot(hist_t/mids**2, mids, c='b', label='T')
                plt.plot(hist_tot/mids**2, mids, c='k')
                plt.xlabel('$g(\omega)/\omega^2$', fontsize='x-large')
                plt.xlim(xlim_bp)
            plt.ylabel('Frequency (THz)', fontsize='x-large')
            plt.ylim(ylim)
            plt.subplots_adjust(left=0.40, bottom=0.25, right=0.60, top=0.752, wspace=0.2, hspace=0.2)
            plt.tick_params(axis="both",direction="in", labelsize='x-large')
            ax.xaxis.set_major_locator(plt.MaxNLocator(1))
            plt.legend(handlelength=1, fontsize='medium').set_draggable(True)
            plt.grid(alpha=0.5)
        plt.show()
