#!/usr/bin/env python
import numpy as np
from subprocess import call
from ase.io import read, write

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
        self.atoms = read(atoms_file, -1)
        self.dynmat = None
        self.w2 = None
        self.w_plot = None
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

    def get_A_mat(self):
        r"""

                                    S_{ij}
        A_{i,j\alpha} (i != j) = ------------ * (r_{j\alpha} - c_{ij\alpha})
                                  V_i r_{ij}

        A_{i,j\alpha} (i == j) = - (sum)_{j != i} A_{i,j\alpha}

        where c_{ij\alpha} is the centroid of the facet polygon S_{ij}.

        """
        if self.A is None:
            npy_name = 'A_mat-{}.npy'.format(self.atoms_file)
            try:
                A = np.load(npy_name)
            except:
                do_calc = True
                print('A matrix file NOT loaded.'.format(npy_name), flush=True)
            else:
                do_calc = False
                print('A matrix file "{}" loaded.'.format(npy_name), flush=True)

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

    def _get_w_plot(
        self,
        nan_to_zero,
        imag_to_minus,
        ):
        if self.w_plot is None:
            from vibration_solver import get_w_from_w2
            self.w_plot = get_w_from_w2(
                self.w2,
                nan_to_zero=nan_to_zero,
                imag_to_minus=imag_to_minus,
                )

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
        self._get_w_plot(nan_to_zero, plot_imaginary)

        hist_l, edges = np.histogram(self.w_plot, bins=nbins, weights=np.linalg.norm(self.ul, axis=-1)**2)
        hist_t, edges = np.histogram(self.w_plot, bins=nbins, weights=np.linalg.norm(self.ut, axis=-1)**2)
        hist_tot, edges = np.histogram(self.w_plot, bins=nbins, density=True)
        mids = (edges[0:-1] + edges[1:]) /2.
        dw = mids[1] -mids[0]
        hist_l /= len(self.w_plot) *dw
        hist_t /= len(self.w_plot) *dw

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

    def get_atomic_VDOS_and_AM(
        self,
        freq_range=(-np.inf, np.inf),
        reduced=True,
        show_2d=True,
        ):
        """
        freq_range (tuple of the form (float1, float2))
            - Frequency range to get VDOS and AM projection.
            - Must be (float1 < float2)
            - In THz unit
        reduced (bool)
            - Reduced VDOS (=VDOS /w^2).
        show_2d (bool)
            - 2-dimensionized images will be written

        Return
            - VDOS will be saved as vdos, vdos_l, vdos_t and AAM as aam in atoms._calc.results dict
            - Real part of sum of eigenvertors will be saved as forces.
        """

        fname = 'vdos-from{}to{}THz-{}.traj'.format(*freq_range, self.atoms_file)
        if reduced:
            fname = 'r' + fname

        try:
            alist = read(fname, ':')
        except:
            do_calc = True
            print('VDOS and AAM file NOT loaded.'.format(fname), flush=True)
        else:
            do_calc = False
            print('VDOS and AAM file "{}" loaded.'.format(fname), flush=True)

        if do_calc:
            # VDOS part
            self.LT_decompose()
            self._get_w_plot(nan_to_zero=True, imag_to_minus=True)

            mask = freq_range[0] < self.w_plot
            mask *= self.w_plot < freq_range[1]
            n_mode = np.sum(mask)

            if reduced:
                wei = np.expand_dims(self.w_plot[mask]**2, axis=1)
            else:
                wei = np.ones((n_mode, 1))
            vdos = np.sum(np.linalg.norm(self.eps[mask].reshape(n_mode, len(self.atoms), 3), axis=-1)**2 /wei, axis=0)
            vdos_l = np.sum(np.linalg.norm(self.ul[mask].reshape(n_mode, len(self.atoms), 3), axis=-1)**2 /wei, axis=0)
            vdos_t = np.sum(np.linalg.norm(self.ut[mask].reshape(n_mode, len(self.atoms), 3), axis=-1)**2 /wei, axis=0)

            # AAM part
            from pam import AAM
            # aam.shape == (n_mode, len(atoms), 3)
            aam = AAM(np.expand_dims(self.eps[mask], axis=0))[0]
            aam = np.sum(aam, axis=0)
            # print(aam.tolist())

            # mode_sum part
            mode_sum = np.sum(self.eps[mask], axis=0).reshape(len(self.atoms), 3)

            alist = []
            atoms = self.atoms.copy()
            from ase.calculators.singlepoint import SinglePointCalculator
            atoms._calc = SinglePointCalculator(atoms, vdos=vdos, vdos_l=vdos_l, vdos_t=vdos_t, aam=aam, forces=mode_sum)
            alist = [atoms]

            if show_2d:
                pos = self.atoms.get_positions()
                from atom_intensity import get_atom_inten_hist, get_peaks
                direcs = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
                for i in range(3):
                    i_hist, edges = get_atom_inten_hist(
                        [self.atoms],
                        direcs[i],
                        nbins=200,
                        )
                    peaks, _ = get_peaks(
                        np.sum(i_hist, axis=0),
                        (edges[0:-1] + edges[1:]) /2.,
                        0.2,
                        )
                    peak_edges = (peaks[0:-1] + peaks[1:])/2.
                    cls = np.digitize(pos[:, i], peak_edges)
                    for j in range(len(peaks)):
                        atoms = self.atoms.copy()[cls == j]
                        atoms._calc = SinglePointCalculator(
                            atoms,
                            vdos=vdos[cls == j],
                            vdos_l=vdos_l[cls == j],
                            vdos_t=vdos_t[cls == j],
                            aam=aam[cls == j],
                            forces=mode_sum[cls == j],
                            )
                        alist.append(atoms)

            write(fname, alist)

        from ase.visualize import view
        view(alist)

        # Histogram
        aam_norm = np.linalg.norm(alist[0].calc.results['aam'], axis=-1)
        hist, bins = np.histogram(aam_norm, 1000)
        mids = (bins[1:] + bins[:-1]) /2.
        from matplotlib import pyplot as plt
        plt.plot(mids, hist, c='k')
        plt.xlabel('Atomic angular momentum (eV ps)', fontsize='x-large')
        plt.ylabel('Population', fontsize='x-large')
        plt.title('mean: {:.4e}, std: {:.4e}'.format(np.mean(aam_norm), np.std(aam_norm)), fontsize='large')
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.subplots_adjust(left=0.20, bottom=0.20, right=0.80, top=0.80, wspace=0.2, hspace=0.2)
        plt.grid(alpha=0.5)

        # Histogram (log)
        aam_norm_log = np.log(aam_norm)
        hist_log, bins_log = np.histogram(aam_norm_log, 50)
        mids_log = (bins_log[1:] + bins_log[:-1]) /2.
        plt.figure()
        plt.plot(np.exp(mids_log), hist_log, c='k')
        plt.xscale('log')
        plt.xlabel('Atomic angular momentum (eV ps)', fontsize='x-large')
        plt.ylabel('Population', fontsize='x-large')
        plt.title('mean: {:.4e}, std: {:.4e}'.format(np.mean(aam_norm), np.std(aam_norm)), fontsize='large')
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.subplots_adjust(left=0.20, bottom=0.20, right=0.80, top=0.80, wspace=0.2, hspace=0.2)
        plt.grid(alpha=0.5)
        plt.show()
