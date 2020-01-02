# FOPTICS
Calculate optical response from multilayer thin film. 
This is a `python3` implementation similar to **Passler2017** and **Vorwerk2018** (see [reference](#reference)).

The semi-infinite thin film is described by a Carteasian coordinate where
 the origin sits at the medium/film surface defining xy plane. The unit vector
  pointing to substrate from the medium is +**z**. The incident plane is xz
  and the incident light is defined by theta, the angle between its wavevector 
  **k** and +**z**. The polarization of incident light is defined by sigma,
  the angle between its projection on xy plane and +**x**.
 
Assuming the medium is isotropic, the inputs for calculating optical response 
from aformentioned system are:
- **theta**: in radians, the incident angle
- **sigma**: in radians, the polarization angle
- **thickness_list**: in nm, the thickness for each layer
- **eps_medium**: real number, the relative dielectric constant of isotropic medium
- **w_list**: a list of frequencies in eV
- **eps_list**: a list of 3x3 matrices, the dielectric tensors for each layer at certain **w**
- **mu**: real number, the magnetic permeability of isotropic medium.

Additionally,
- **Euler_alpha**, **Euler_beta**, **Euler_gamma**: in radians, the Euler angles used to convert **eps_list**
to current Cartesians. 

Usage
---
1. you need `numpy`.
2. prepare input epsilon files:

    The code reads **eps_list** and **w_list** from 2 table-like files containing 
    either the imaginary or the real part of epsilon. 
    Both of them have the format:
    ```bash
    # w  xx  xy  xz  yx  yy  yz  zx  zy  zz
    ```
    so one needs to make sure the 2 inputs share the identical 1st column.
    
    For VASP users, there is a script `vasp_eps.sh` here can be used to extract the 2 input files
    from `vasprun.xml`.
3. see `test.py` for usage

Issues
---
As pointed out in the erratum of **Passler2017**, total transmittance cannot be calculated directly from the transmission coefficients unless the substrate is vacuum. It seems this restriction is ignored in **Vorwerk2018**.


Reference 
---
**Berreman1971**: Berreman, Dwight W. "Optics in stratified and anisotropic media: 4× 4-matrix formulation." Josa 62.4 (1972): 502-510.

**Yeh1979**: Yeh, Pochi. "Optics of anisotropic layered media: a new 4× 4 matrix algebra." Surface Science 96.1-3 (1980): 41-53.

**Xu2000**: Xu, W., L. T. Wood, and T. D. Golding. "Optical degeneracies in anisotropic layered media: treatment of singularities in a 4× 4 matrix formalism." Physical Review B 61.3 (2000): 1740.

**Passler2017**: Passler, Nikolai Christian, and Alexander Paarmann. "Generalized 4× 4 matrix formalism for light propagation in anisotropic stratified media: study of surface phonon polaritons in polar dielectric heterostructures." JOSA B 34.10 (2017): 2128-2139.

**Vorwerk2018**: Vorwerk, Christian, Caterina Cocchi, and Claudia Draxl. "LayerOptics: Microscopic modeling of optical coefficients in layered materials." Computer Physics Communications 201 (2016): 119-125.
