module constants
  double precision, parameter :: conv_ev = 27.2113860217
  double precision, parameter :: conv_bohr = 0.52917721090380
  double precision, parameter :: K_PARAMETER = 1.75

  ! Atom name
  character(len=2), parameter::atomname(110) = (/'H ', 'He', 'Li', 'Be', 'B ', 'C ', &
                                                 'N ', 'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', &
                                                 'Cl', 'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', &
                                                 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
                                                 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', &
                                                 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', &
                                                 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', &
                                                 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', &
                                                 'Ir', 'Pt', 'Au', 'Hg', 'Ti', 'Pb', 'Bi', 'Po', 'At', 'Rn', &
                                                 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', &
                                                 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', &
                                                 'Bh', 'Hs', 'Mt', 'Ds'/)

  ! Atomic radii were obtained from
  ! https://doi.org/10.1002/chem.201602949
  ! as compiled by the mendeleev, and converted to a.u.
  ! see script radii.py
  double precision, parameter :: ATOMIC_RADII(110) = (/ &
                                                     0.815d0, 0.709d0, 1.164d0, 1.159d0, 1.085d0, 1.005d0, 0.947d0, 0.905d0, &
                                                     0.863d0, 0.826d0, 1.191d0, 1.270d0, 1.265d0, 1.228d0, 1.180d0, 1.132d0, &
                                                     1.090d0, 1.042d0, 1.238d0, 1.429d0, 1.392d0, 1.360d0, 1.334d0, 1.233d0, &
                                                     1.281d0, 1.254d0, 1.233d0, 1.212d0, 1.148d0, 1.175d0, 1.233d0, 1.238d0, &
                                                     1.222d0, 1.185d0, 1.159d0, 1.122d0, 1.270d0, 1.476d0, 1.450d0, 1.423d0, &
                                                     1.328d0, 1.291d0, 1.334d0, 1.254d0, 1.233d0, 1.138d0, 1.191d0, 1.259d0, &
                                                     1.302d0, 1.312d0, 1.302d0, 1.281d0, 1.259d0, 1.228d0, 1.318d0, 1.550d0, &
                                                     1.503d0, 1.492d0, 1.513d0, 1.503d0, 1.498d0, 1.482d0, 1.482d0, 1.466d0, &
                                                     1.461d0, 1.455d0, 1.445d0, 1.439d0, 1.434d0, 1.466d0, 1.429d0, 1.397d0, &
                                                     1.365d0, 1.339d0, 1.318d0, 1.291d0, 1.270d0, 1.217d0, 1.196d0, 1.212d0, &
                                                     1.281d0, 1.318d0, 1.323d0, 1.323d0, 1.307d0, 1.286d0, 1.365d0, 1.545d0, &
                                                     1.550d0, 1.524d0, 1.508d0, 1.498d0, 1.487d0, 1.471d0, 1.461d0, 1.461d0, &
                                                     1.482d0, 1.482d0, 1.482d0, 1.482d0, 1.482d0, 1.482d0, 1.482d0, 1.482d0, &
                                                     1.482d0, 1.482d0, 1.482d0, 1.482d0, 1.482d0, 1.482d0/)

end module constants
