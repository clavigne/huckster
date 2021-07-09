"""Build atomic data to be included in the bond detection program."""
import pyscf
from pyscf import gto
from pyscf.scf import atom_hf
import numpy as np
import shutil
import os.path

# Elements and basis fuctions
atoms = dict([(index, 'sto-3g') for index in range(1,54)])
MAX_NATOMS = 1000

# parameters
ptr = 0

print("Welcome to the atomic SCF dumper --- powered by pyscf.")
bas = []
env = []

aos = np.zeros((110,5), dtype=np.int)
at_off = 1

params = ""

for key, val in atoms.items():
    print("===================================================================================")
    print(f"Atom:  {key}")
    print(f"Basis: {val}")
    print(f"env -> {ptr}")
    mol = gto.Mole()
    mol.verbose = 0
    mol.output = None
    mol.atom =[[key, (0.0, 0.0, 0.0)]]
    mol.basis = {key: val}
    mol.spin = None
    mol.build()

    print(f"starting atomic scf...")
    atm_scf = atom_hf.get_atm_nrhf(mol)
    print("done!")

    # Save all the data related to the atomic orbitals
    e_tot, mo_energy, mo_coeff, mo_occ = list(atm_scf.values()).pop()
    print("    AO energies:", mo_energy)
    print("AO coefficients:", mo_coeff[0])
    for i in range(1, len(mo_energy)):
        print("                ", mo_coeff[i])
    print(" AO occupations:", mo_occ)

    mo_energy_ptr = ptr
    env += [el for el in mo_energy]
    ptr += len(mo_energy)

    mo_coeff_ptr = ptr
    for vec in mo_coeff:
        env += [el for el in vec]
        ptr += len(vec)

    mo_occ_ptr = ptr
    env += [el for el in mo_occ]
    ptr += len(mo_occ)

    aos[key-1,:] = np.array([at_off, len(mo_energy), mo_energy_ptr, mo_coeff_ptr, mo_occ_ptr])

    # Only one element so that's easy
    atombasis = list(mol._basis.values()).pop()

    # this next part is basically verbatim from
    # https://github.com/pyscf/pyscf/blob/53e2069b4a3a2e0616bdf4d8c2e3f898c10a8330/pyscf/gto/mole.py#L876
    print("converting to libcint formats...")
    b, e = gto.make_bas_env(atombasis, atom_id=key, ptr=ptr)
    print("done!")

    # Now we have "prototype" bas and env arrays to give to the fortran code.
    # Note that we replaced the atom index in the basis vector by the atom's Z
    # value, so that we can easily find them.
    bas += [b[:]]

    env += list(e)
    ptr += len(e)
    at_off += len(b)

# Done
print("===================================================================================")
bas = np.array(np.concatenate(bas).T)
aos = np.array(aos).T
env = np.array(env)
print('Size of bas array:', bas.shape)
print('Size of aos array:', aos.shape)
print('Size of env array:', env.shape)

def make_array_2d(name, data, fstr):
    out = f'{name} = reshape((/ '
    k = 0
    n = np.product(data.shape)
    while True:
        for i in range(10):
            out += fstr % (data.T.flat[k]) 
            k += 1
            if k == n:
                break
            else:
                out += ','

        if k == n:
            break
        else:
            out += '&\n    '
    return out + f' /), (/ {data.shape[0]}, {data.shape[1]} /) )\n'

def make_array_1d(name, data, fstr):
    out = f'{name} = (/ '
    k = 0
    n = np.product(data.shape)
    while True:
        for i in range(8):
            out += fstr % (data.flat[k]) 
            k += 1
            if k == n:
                break
            else:
                out += ','

        if k == n:
            break
        else:
            out += '&\n    '
    return out + f' /)\n'
        


from preamble import preamble
print("dumping data in fortran files.")
sizeof_env = env.shape[0] + 3 * MAX_NATOMS  # we store the coordinates here too
with open('data-header.f90', 'wb') as f:
    f.write(preamble.encode('ascii'))
    f.write(f'integer, parameter :: SIZEOF_BAS = {bas.shape[1]} \n'.encode('ascii'))
    f.write(f'integer, parameter :: SIZEOF_AOS = {aos.shape[1]} \n'.encode('ascii'))
    f.write(f'integer, parameter :: SIZEOF_ENV = {sizeof_env} \n'.encode('ascii'))
    f.write(f'integer, parameter :: ENV_COORD_PTR = {env.shape[0]} \n'.encode('ascii'))
    f.write(f'integer, parameter :: MAX_NATOMS = {MAX_NATOMS} \n'.encode('ascii'))

    f.write(f'integer :: pbas(BAS_SLOTS, SIZEOF_BAS)\n'.encode('ascii'))
    f.write(f'integer :: paos(AOS_SLOTS, SIZEOF_AOS)\n'.encode('ascii'))
    f.write(f'double precision :: ENV(SIZEOF_ENV)\n'.encode('ascii'))


with open('data-content.f90', 'wb') as f:
    f.write("""
! DATA FROM ATOMIC SCF CALCULATIONS
""".encode('ascii'))
    f.write(make_array_2d('pbas', bas, '%i').encode('ascii'))
    f.write(make_array_2d('paos', aos, '%i').encode('ascii'))

    # env is so big that we need to break it up
    blocks = 2048
    ll = 1
    ul = blocks
    while True:
        f.write(make_array_1d(f'env({ll}:{ul-1})', env[ll:ul], '%f').encode('ascii'))
        if ul == env.shape[0]:
            break

        ll += blocks
        ul += blocks

        if ul > env.shape[0]:
            ul = env.shape[0]

# Copy to destination
shutil.copy('data-header.f90', os.path.join('..', 'include', 'data-header.f90'))
shutil.copy('data-content.f90', os.path.join('..', 'include', 'data-content.f90'))

# Finally save the basis for tests.
import pickle
with open('basis.pkl', 'wb') as f:
    pickle.dump(atoms, f)

