import numpy as np
import math
import sys

"""
reference:
**Berreman1971**: Berreman, Dwight W. "Optics in stratified and anisotropic media: 4× 4-matrix formulation." Josa 62.4 (1972): 502-510.

**Yeh1979**: Yeh, Pochi. "Optics of anisotropic layered media: a new 4× 4 matrix algebra." Surface Science 96.1-3 (1980): 41-53.

**Xu2000**: Xu, W., L. T. Wood, and T. D. Golding. "Optical degeneracies in anisotropic layered media: treatment of singularities in a 4× 4 matrix formalism." Physical Review B 61.3 (2000): 1740.

**Passler2017**: Passler, Nikolai Christian, and Alexander Paarmann. "Generalized 4× 4 matrix formalism for light propagation in anisotropic stratified media: study of surface phonon polaritons in polar dielectric heterostructures." JOSA B 34.10 (2017): 2128-2139.

**Vorwerk2018**: Vorwerk, Christian, Caterina Cocchi, and Claudia Draxl. "LayerOptics: Microscopic modeling of optical coefficients in layered materials." Computer Physics Communications 201 (2016): 119-125.
"""


class Layer_Calculator:

    def __init__(self, w, theta, eps_medium, mu, thickness):
        self.w = w
        self.zeta = np.sqrt(eps_medium) * np.sin(theta)  # conserve throughout
        self.mu = mu
        self.theta = theta
        self.thickness = thickness * 1e-9

    @staticmethod
    def eq13(v):
        """
        this is Passler2017_eq13
        :param v: eigenmode psi in Cartesian
        :return:
        """
        if (abs(v[0]) ** 2 + abs(v[2]) ** 2) == 0:
            return 0.0
        return abs(v[0]) ** 2 / (abs(v[0]) ** 2 + abs(v[2]) ** 2)

    def eq11(self, eps):
        """
        solve Passler2017_eq11
        :param eps:
        :return:
        """
        mu = self.mu
        zeta = self.zeta
        M = np.zeros((6, 6), dtype='complex')
        M[0:3, 0:3] = np.array(eps)
        M[3:6, 3:6] = np.eye(3) * mu
        d = M[2, 2] * M[5, 5] - M[2, 5] * M[5, 2]
        # a matrix, Berreman1971_eq21
        a = np.zeros((6, 6), dtype='complex')
        a[2, 0] = (M[5, 0] * M[2, 5] - M[2, 0] * M[5, 5]) / d
        a[2, 1] = ((M[5, 1] - zeta) * M[2, 5] - M[2, 1] * M[5, 5]) / d
        a[2, 3] = (M[5, 3] * M[2, 5] - M[2, 3] * M[5, 5]) / d
        a[2, 4] = (M[5, 4] * M[2, 5] - (M[2, 4] + zeta) * M[5, 5]) / d
        a[5, 0] = (M[5, 2] * M[2, 0] - M[2, 2] * M[5, 0]) / d
        a[5, 1] = (M[5, 2] * M[2, 1] - M[2, 2] * (M[5, 1] - zeta)) / d
        a[5, 3] = (M[5, 2] * M[2, 3] - M[2, 2] * M[5, 3]) / d
        a[5, 4] = (M[5, 2] * (M[2, 4] + zeta) - M[2, 2] * M[5, 4]) / d
        # S matrix, Berreman1971_eq24
        S = np.zeros((4, 4), dtype='complex')
        S[0, 0] = M[0, 0] + M[0, 2] * a[2, 0] + M[0, 5] * a[5, 0]
        S[0, 1] = M[0, 1] + M[0, 2] * a[2, 1] + M[0, 5] * a[5, 1]
        S[0, 2] = M[0, 3] + M[0, 2] * a[2, 3] + M[0, 5] * a[5, 3]
        S[0, 3] = M[0, 4] + M[0, 2] * a[2, 4] + M[0, 5] * a[5, 4]
        S[1, 0] = M[1, 0] + M[1, 2] * a[2, 0] + (M[1, 5] - zeta) * a[5, 0]
        S[1, 1] = M[1, 1] + M[1, 2] * a[2, 1] + (M[1, 5] - zeta) * a[5, 1]
        S[1, 2] = M[1, 3] + M[1, 2] * a[2, 3] + (M[1, 5] - zeta) * a[5, 3]
        S[1, 3] = M[1, 4] + M[1, 2] * a[2, 4] + (M[1, 5] - zeta) * a[5, 4]
        S[2, 0] = M[3, 0] + M[3, 2] * a[2, 0] + M[3, 5] * a[5, 0]
        S[2, 1] = M[3, 1] + M[3, 2] * a[2, 1] + M[3, 5] * a[5, 1]
        S[2, 2] = M[3, 3] + M[3, 2] * a[2, 3] + M[3, 5] * a[5, 3]
        S[2, 3] = M[3, 4] + M[3, 2] * a[2, 4] + M[3, 5] * a[5, 4]
        S[3, 0] = M[4, 0] + (M[4, 2] + zeta) * a[2, 0] + M[3, 5] * a[5, 0]
        S[3, 1] = M[4, 1] + (M[4, 2] + zeta) * a[2, 1] + M[3, 5] * a[5, 1]
        S[3, 2] = M[4, 3] + (M[4, 2] + zeta) * a[2, 3] + M[3, 5] * a[5, 3]
        S[3, 3] = M[4, 4] + (M[4, 2] + zeta) * a[2, 4] + M[3, 5] * a[5, 4]
        # Delta matrix, Berreman1971_eq24
        Dt = np.zeros((4, 4), dtype='complex')
        Dt[0, 0] = S[3, 0]
        Dt[0, 1] = S[3, 3]
        Dt[0, 2] = S[3, 1]
        Dt[0, 3] = - S[3, 2]
        Dt[1, 0] = S[0, 0]
        Dt[1, 1] = S[0, 3]
        Dt[1, 2] = S[0, 1]
        Dt[1, 3] = - S[0, 2]
        Dt[2, 0] = -S[2, 0]
        Dt[2, 1] = -S[2, 3]
        Dt[2, 2] = -S[2, 1]
        Dt[2, 3] = S[2, 2]
        Dt[3, 0] = S[1, 0]
        Dt[3, 1] = S[1, 3]
        Dt[3, 2] = S[1, 1]
        Dt[3, 3] = - S[1, 2]
        qs, psis = np.linalg.eig(Dt)
        return qs, psis, a

    def eq14(self, qs, psis, mat_a):
        """
        use Passler2017_eq14 to sort eigenmodes and calculate corresponding Poynting vector
        :param qs:
        :param psis:
        :param mat_a:
        :return:
        """

        # first use eq12 to classify
        transmode = []
        reflmode = []
        if np.any(np.iscomplex(qs)):
            zeroidx = []
            for i in range(0, 4):
                qval = qs[i]
                if qval.imag > 0:
                    transmode.append(i)
                elif qval.imag < 0:
                    reflmode.append(i)
                else:
                    zeroidx.append(i)

        else:
            zeroidx = []
            for i in range(0, 4):
                qval = qs[i]
                if qval.real > 0:
                    transmode.append(i)
                elif qval.real < 0:
                    reflmode.append(i)
                else:
                    zeroidx.append(i)

        if len(transmode) == 2 and len(reflmode) == 2:
            pass
        elif len(transmode) == 1 and len(reflmode) == 1:
            transmode += [zeroidx[0]]
            reflmode += [zeroidx[1]]
        elif len(transmode) == 0 and len(reflmode) == 0 and len(zeroidx) == 4:
            transmode = zeroidx[:2]
            reflmode = zeroidx[2:]
        else:
            sys.exit('something is wrong at w={} with qs: {}'.format(self.w, qs))

        Py = np.zeros((3, 4), dtype='complex')
        # calculate Poynting vector for each Eigenmode psi, Py[:, i] is the vector
        for i in range(0, 4):
            Ex = psis[0, i]
            Ey = psis[2, i]
            Hx = -psis[3, i]
            Hy = psis[1, i]

            Ez = mat_a[2, 0] * Ex + mat_a[2, 1] * Ey + mat_a[2, 3] * Hx + mat_a[2, 4] * Hy
            Hz = mat_a[5, 0] * Ex + mat_a[5, 1] * Ey + mat_a[5, 3] * Hx + mat_a[5, 4] * Hy
            Pyx = Ey * Hz - Ez * Hy
            Pyy = Ez * Hx - Ex * Hz
            Pyz = Ex * Hy - Ey * Hx
            Py[0, i] = Pyx
            Py[1, i] = Pyy
            Py[2, i] = Pyz

        Cp1 = self.eq13(Py[:, transmode[0]])
        Cp2 = self.eq13(Py[:, transmode[1]])

        if abs(Cp1 - Cp2) > 1e-15:  # birefringence -> use Poynting vector
            if Cp2 > Cp1:
                transmode[0], transmode[1] = transmode[1], transmode[0]
            if self.eq13(Py[:, reflmode[1]]) > self.eq13(Py[:, reflmode[0]]):
                reflmode[0], reflmode[1] = reflmode[1], reflmode[0]
        else:  # no birefringence -> use s-pol/p-pol
            if self.eq13(psis[:, transmode[1]]) > self.eq13(psis[:, transmode[0]]):
                transmode[0], transmode[1] = transmode[1], transmode[0]
            if self.eq13(psis[:, reflmode[1]]) > self.eq13(psis[:, reflmode[0]]):
                reflmode[0], reflmode[1] = reflmode[1], reflmode[0]

        psis = [
            psis[:, transmode[0]],
            psis[:, transmode[1]],
            psis[:, reflmode[0]],
            psis[:, reflmode[1]],
        ]

        qs = [qs[transmode[0]],
              qs[transmode[1]],
              qs[reflmode[0]],
              qs[reflmode[1]]]

        Py = [Py[:, transmode[0]],
              Py[:, transmode[1]],
              Py[:, reflmode[0]],
              Py[:, reflmode[1]]]

        return psis, qs, Py

    def eq20(self, qs, eps):
        """
        use Passler2017_eq20, in the matlab code of Passler2017 there seems to be a typo for gamma13
        the implementation here observes the paper and also Xu2000
        :param qs:
        :param eps:
        :return: gamma matrix, representing polarization vectors
        """
        zeta = self.zeta
        mu = self.mu
        eps = mu * eps
        gamma = np.zeros((4, 3), dtype='complex')
        gamma[0, 0] = 1
        gamma[1, 1] = 1
        gamma[3, 1] = 1
        gamma[2, 0] = -1
        # often needed denominator
        den_mueps33mq02 = (eps[2, 2] - zeta ** 2)

        if abs(qs[0] - qs[1]) < 1e-15:
            gamma12 = 0.0
            gamma13 = - (eps[2, 0] + zeta * qs[0]) / den_mueps33mq02
            gamma21 = 0.0
            gamma23 = - eps[2, 1] / den_mueps33mq02
        else:
            gamma12_de = den_mueps33mq02 * (eps[1, 1] - zeta ** 2 - qs[0] ** 2) - eps[1, 2] * eps[2, 1]
            if gamma12_de:
                gamma12 = (eps[1, 2] * (eps[2, 0] + zeta * qs[0]) - eps[1, 0] * den_mueps33mq02) / (
                        den_mueps33mq02 * (eps[1, 1] - zeta ** 2 - qs[0] ** 2) - eps[1, 2] * eps[2, 1])
            else:
                gamma12 = 0.0

            gamma13_de = den_mueps33mq02 * gamma12
            if gamma13_de:
                gamma13 = - (eps[2, 0] + zeta * qs[0]) / den_mueps33mq02 - eps[2, 1] / den_mueps33mq02 * gamma12
                # gamma13 = - (eps[2, 0] + zeta * qs[0]) / den_mueps33mq02 + eps[2, 1] / den_mueps33mq02 * gamma12
            else:
                gamma13 = - (eps[2, 0] + zeta * qs[0]) / den_mueps33mq02

            gamma21_de = den_mueps33mq02 * (eps[0, 0] - qs[1] ** 2) - (eps[0, 2] + zeta * qs[1]) * (
                    eps[2, 0] + zeta * qs[1])
            if gamma21_de:
                gamma21 = (eps[2, 1] * (eps[0, 2] + zeta * qs[1]) - eps[0, 1] * den_mueps33mq02) / (
                        den_mueps33mq02 * (eps[0, 0] - qs[1] ** 2) - (eps[0, 2] + zeta * qs[1]) * (
                        eps[2, 0] + zeta * qs[1]))
            else:
                gamma21 = 0.0

            gamma23_de = den_mueps33mq02 * gamma21
            if gamma23_de:
                gamma23 = - (eps[2, 0] + zeta * qs[1]) / den_mueps33mq02 * gamma21 - eps[2, 1] / den_mueps33mq02
            else:
                gamma23 = - eps[2, 1] / den_mueps33mq02

        if abs(qs[2] - qs[3]) < 1e-15:
            gamma32 = 0.0
            gamma33 = (eps[2, 0] + zeta * qs[2]) / den_mueps33mq02
            gamma41 = 0.0
            gamma43 = - eps[2, 1] / den_mueps33mq02
        else:
            if den_mueps33mq02 * (eps[1, 1] - zeta ** 2 - qs[2] ** 2) - eps[1, 2] * eps[2, 1]:
                gamma32 = (eps[1, 0] * den_mueps33mq02 - eps[1, 2] * (eps[2, 0] + zeta * qs[2])) / (
                        den_mueps33mq02 * (eps[1, 1] - zeta ** 2 - qs[2] ** 2) - eps[1, 2] * eps[2, 1])
            else:
                gamma32 = 0.0

            if den_mueps33mq02 * gamma32:
                gamma33 = (eps[2, 0] + zeta * qs[2]) / den_mueps33mq02 + eps[2, 1] / den_mueps33mq02 * gamma32
            else:
                gamma33 = (eps[2, 0] + zeta * qs[2]) / den_mueps33mq02

            if den_mueps33mq02 * (eps[0, 0] - qs[3] ** 2) - (eps[0, 2] + zeta * qs[3]) * (eps[2, 0] + zeta * qs[3]):
                gamma41 = (eps[2, 1] * (eps[0, 2] + zeta * qs[3]) - eps[0, 1] * den_mueps33mq02) / (
                        den_mueps33mq02 * (eps[0, 0] - qs[3] ** 2) - (eps[0, 2] + zeta * qs[3]) * (
                        eps[2, 0] + zeta * qs[3]))
            else:
                gamma41 = 0.0

            if den_mueps33mq02 * gamma41:
                gamma43 = - (eps[2, 0] + zeta * qs[3]) / den_mueps33mq02 * gamma41 - eps[2, 1] / den_mueps33mq02
            else:
                gamma43 = - eps[2, 1] / den_mueps33mq02

        gamma[0, 1] = gamma12
        gamma[0, 2] = gamma13
        gamma[1, 0] = gamma21
        gamma[1, 2] = gamma23
        gamma[2, 1] = gamma32
        gamma[2, 2] = gamma33
        gamma[3, 0] = gamma41
        gamma[3, 2] = gamma43
        return gamma

    def eq26(self, gamma, qs, ):
        """
        calculate transfer matrix for layer as Passler2017_eq26, notice this has not been transformed by Passler2017_eq32
        :param gamma:
        :param qs:
        :return:
        """
        Ai = np.zeros((4, 4), dtype='complex')
        Ki = np.zeros((4, 4), dtype='complex')
        Ai[0] = gamma[:, 0]
        Ai[1] = gamma[:, 1]
        for k in range(4):
            Ai[2, k] = 1 / self.mu * qs[k] * gamma[k, 0] - self.zeta * gamma[k, 2]
            Ai[3, k] = 1 / self.mu * qs[k] * gamma[k, 1]
            expidx = np.complex(-1j * 2 * np.pi * self.w * 8065.54429 * 1e2 * qs[k] * self.thickness)
            Ki[k, k] = np.exp(expidx)
        Ti = Ai @ Ki @ np.linalg.inv(Ai)
        return Ai, Ki, Ti


class Calculator:

    def __init__(self, w, theta, eps_medium, mu, thickness_list, sigma, eps_list, Euler_alpha=0.0, Euler_beta=0.0,
                 Euler_gamma=0.0):
        """
        init with a certain w, eps_list is for layers, not for a list of ws
        :param w: in eV
        :param theta: in radians
        :param eps_medium: relative dielectric
        :param mu: relative permeability
        :param thickness_list: in nm
        :param sigma: in radians
        :param eps_list: relative dielectric for layers, not for a list of ws, eps_list[i] is the eps matrix of the i-1th layer at a certain w
        :param Euler_alpha: in radians
        :param Euler_beta: in radians
        :param Euler_gamma: in radians
        """
        self.w = w
        # self.zeta = np.sqrt(eps_medium) * np.sin(theta)  # conserve throughout
        self.mu = mu
        self.theta = theta
        self.thickness_list = thickness_list
        self.sigma = sigma

        rotmat = self.euler_rotmat(Euler_alpha, Euler_beta, Euler_gamma)
        rotmat_inv = np.linalg.inv(rotmat)
        for i in range(len(eps_list)):
            eps_list[i] = rotmat @ eps_list[i] @ rotmat_inv
        self.eps_list = eps_list

        self.eps_medium = eps_medium
        self.theta = theta

    @staticmethod
    def euler_rotmat(alpha, beta, gamma):
        """
        euler rotation matrix, all in radians
        :param alpha: rot x
        :param beta: rot y
        :param gamma: rot z
        :return:
        """
        c1 = math.cos(alpha)
        c2 = math.cos(beta)
        c3 = math.cos(gamma)
        s1 = math.sin(alpha)
        s2 = math.sin(beta)
        s3 = math.sin(gamma)
        rotmat = np.zeros((3, 3))
        rotmat[0, 0] = c1 * c3 - c2 * s2 * s3
        rotmat[0, 1] = -c1 * s3 - c2 * c3 * s1
        rotmat[0, 2] = s1 * s2
        rotmat[1, 0] = c3 * s1 + c1 * c2 * s3
        rotmat[1, 1] = c1 * c2 * c3 - s1 * s3
        rotmat[1, 2] = -c1 * s2
        rotmat[2, 0] = s2 * s3
        rotmat[2, 1] = c3 * s2
        rotmat[2, 2] = c2
        return rotmat

    def cal(self):
        params = []

        for ilayer in range(len(self.thickness_list)):
            thickness = self.thickness_list[ilayer]
            eps = self.eps_list[ilayer]

            lc = Layer_Calculator(self.w, self.theta, self.eps_medium, self.mu, thickness)
            qs, psis, mat_a = lc.eq11(eps)
            psis, qs, Py = lc.eq14(qs, psis, mat_a)
            gamma = lc.eq20(qs, eps)
            Ai, Ki, Ti = lc.eq26(gamma, qs)
            params.append([Ai, Ki, Ti, Py, gamma])

        A0, K0, T0, Pyi0, gamma0 = params[0]
        Af, Kf, Tf, Pyif, gammaf = params[-1]
        SwapColumns1324 = [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
        SwapColumns1324 = np.array(SwapColumns1324)
        T = np.eye(4)

        for i in range(len(self.eps_list) - 2, 0, -1):
            AI, KI, TI, Pyord, gamma = params[i]
            T = TI @ T
        T = np.linalg.inv(A0) @ T @ Af
        T = SwapColumns1324 @ T @ SwapColumns1324

        # Yeh1979_eq 22, this can only be used after swap
        M11, M12, M13, M14 = T[0]
        M21, M22, M23, M24 = T[1]
        M31, M32, M33, M34 = T[2]
        M41, M42, M43, M44 = T[3]

        As = np.sin(self.sigma)
        Ap = np.cos(self.sigma)
        Cs, Cpol = np.linalg.inv([[M11, M13], [M31, M33]]) @ np.array([As, Ap])
        Bs, Bp = np.linalg.inv([[M21, M23], [M41, M43]]) @ np.array([Cs, Cpol])
        R1 = abs(Bs) ** 2
        R2 = abs(Bp) ** 2
        T1 = abs(Cs) ** 2
        T2 = abs(Cpol) ** 2
        absor = -np.log(T1 + T2)
        return R1, R2, T1, T2, absor

    @staticmethod
    def records_from_file(re_file, im_file):
        re_epsilon = np.loadtxt(re_file, dtype=np.float64)
        im_epsilon = np.loadtxt(im_file, dtype=np.float64)

        records = []
        for i in range(len(re_epsilon)):
            epsilon = dict()
            omega = re_epsilon[i][0]
            elements = re_epsilon[i][1:] + im_epsilon[i][1:] * 1.0j
            epsilon['omega'] = omega
            epsilon['matrix'] = elements.reshape((3, 3))
            records.append(epsilon)
        return records
