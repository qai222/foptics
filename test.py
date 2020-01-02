from foptics import Calculator, np

"""
testing
data for sibse come from vaspwiki, https://cms.mpi.univie.ac.at/wiki/index.php/Dielectric_properties_of_Si_using_BSE
data for bithiophene come from LayerOptics (Vorwerk2018), http://exciting-code.org/carbon-layeroptics-tutorial
"""

sibse_re = './test/sibse.re.epsilon'
sibse_im = './test/sibse.im.epsilon'

bithiophene_re = './test/bithiophene.re.epsilon'
bithiophene_im = './test/bithiophene.im.epsilon'


def test1layer(re, im, theta, eps_medium, eps_substrate, mu, thickness_list, sigma):
    records = Calculator.records_from_file(re, im)
    x = []
    y = []
    for r in records:
        eps = r['matrix']
        w = r['omega']
        eps_list = [eps_medium * np.eye(3, dtype='complex'), eps, eps_substrate * np.eye(3, dtype='complex')]
        cal = Calculator(w, theta, eps_medium, mu, thickness_list, sigma, eps_list)
        R1, R2, T1, T2, absor = cal.cal()
        x.append(w)
        y.append(absor)
    return x, y


import matplotlib.pyplot as plt

theta = 0
sigma = 0
eps_medium = 1
eps_substrate = 1
mu = 1
thickness_list = [0, 0.5, 0]

x1, y1 = test1layer(sibse_re, sibse_im, theta, eps_medium, eps_substrate, mu, thickness_list, sigma)
plt.plot(x1, y1)
plt.title('Si_BSE')
plt.xlabel('w (eV)')
plt.ylabel('abs (a. u.)')
plt.xlim([0, 10])
plt.savefig('./test/sibse.eps')
plt.clf()

for sigma in [0, np.pi / 6, np.pi / 3, np.pi / 2]:
    x2, y2 = test1layer(bithiophene_re, bithiophene_im, theta, eps_medium, eps_substrate, mu, thickness_list, sigma)
    plt.plot(x2, y2)

plt.title('bithiophene')
plt.xlabel('w (eV)')
plt.ylabel('abs (a. u.)')
plt.xlim([3, 5])
plt.savefig('./test/bithiophene.eps')
