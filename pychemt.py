# First attempt at creating a molecule model

# Importando las librerias necesarias
try:
    from LA import *
except:
    raise ImportError, "Couldn't get the necessary libraries."

try:
    import math
    import json
    import array
except:
    raise ImportError, "Interphase: Math and/or JSON could not be found!"

# Abriendo archivo para cargar la tabla periodica a RAM
with open("periodic-table.json") as file_data:
    raw_data = file_data.read()

# Interpretando informacion de la base de datos
data = json.loads(raw_data)

# Creando diccionario de la tabla periodica
PERIODIC_TABLE = {}
for ele in data:
    letras = str(ele["symbol"].upper().split()[0])
    PERIODIC_TABLE[letras] = {}
    llaves = ele.keys()
    for prop in llaves:
        PERIODIC_TABLE[letras][prop] = ele[prop]

# Carga electrica de un electron
K_ELE = 8.987552 * 10 ** -9

# Radio de Bohr
BOHR = 1.0 / 0.5291772083

# Clase atomo para representar cada atomo con sus propias
# coordenadas, simbolo, metodos, etc.
class Atom(object):
    """ Clase para definir bien un atomo y poderlo
    utilizar despues en moleculas. """

    def __init__(self, element="H", x=0.0, y=0.0, z=0.0,
                 charge=0.0, flag=False):
        """ Constructor para crear un atomo nuevo. """
        self.element = element
        self.coords = Matrix([[x, y, z]])
        self.charge = charge
        self.flag = flag

    def __repr__(self):
        """ Metodo para representar el objeto Atom. """
        if self.flag:
            temp = "*"
        else:
            temp = "x"
        text = "<" + self.element + "-" + temp + ">"
        return text

    def set_coordinates(self, x, y, z):
        """ Metodo para actualizar coordenadas simultaneamente. """
        self.coords = Matrix([[x, y, z]])

    def get_coordinates(self):
        """ Metodo para obtener coordenadas en forma de lista. """
        return self.coords[0]


# Clase molecula para representar moleculas y poder calcular
# propiedades de las mismas.
class Molecule(object):
    """ Clase para definir bien una molecula y poderla
    manipular, asi como calcular propiedades de ella. """

    def __init__(self):
        """ Constructor para crear una nueva molecula. """
        self.atoms = []
        self.bonds = []
        self.mol_weight = 0.0
        self.charge = 0.0

    def __repr__(self):
        """ Metodo para representar un objeto molecula. """
        temp = "Atoms:\n\t"
        for a in self.atoms:
            temp += str(a) + "\t"
        return temp

    def add_atoms(self, *atoms):
        """ Metodo para aniadir atomos a la molecula. """
        if len(atoms) == 0:
            raise TypeError, "Molecule: The object added is empty."
        elif len(atoms) == 1:
            if isinstance(atoms[0], list):
                atoms = atoms[0]
        for a in atoms:
            if not isinstance(a, Atom):
                raise TypeError, "Molecule: The object added is not an instance of Atom."
        for a in atoms:
            self.atoms.append(a)
        self.get_mol_weight()
        return True

    def set_bonds(self, *bonds):
        """ Metodo para agregar enlaces a la molecula. """
        if len(bonds) == 0:
            raise TypeError, "Molecule: The object added is empty."
        elif len(bonds) == 1:
            if isinstance(bonds[0], list):
                bonds = bonds[0]
        todos = len(self.atoms)
        self.bonds = [[0] * todos for i in xrange(todos)]
        for b in bonds:
            if not isinstance(b, list) and not isinstance(b, tuple):
                raise TypeError, "Molecule: The given argument is not a list."
            if len(b) != 3:
                raise ValueError, "Molecule: The given argument is incomplete - atom[1], atom[2], number of bonds."
            if (b[0] > todos) or (b[1] > todos) or (b[0] < 0) or \
                    (b[1] < 0) or (not isinstance(b[0], int)) or \
                    (not isinstance(b[1], int)):
                raise ValueError, "Molecule: The atom referenced does not exist."
            self.bonds[b[0]][b[1]] = b[2]
            self.bonds[b[1]][b[0]] = b[2]
        self.bonds = Matrix(self.bonds)
        return True

    def get_mol_weight(self):
        """ Metodo para calcular el peso de la molecula. """
        self.mol_weight = 0.0
        for a in self.atoms:
            symbol = a.element
            self.mol_weight += PERIODIC_TABLE[symbol]["mass"]
        return True

    def get_coords(self):
        """ Metodo para obtener una lista con el id coordenadas
        y carga de la molecula. """
        todos = []
        for a in self.atoms:
            x, y, z = a.get_coordinates()
            todos.append([a.element, x, y, z, a.charge])
        return todos

    def get_center_of_mass(self):
        """ Metodo para obtener el centro de masa de la molecula. """
        atomos = self.get_coords()
        M = self.mol_weight
        centro = Matrix([[0, 0, 0]])
        for a in atomos:
            centro = centro + (Matrix([[a[1], a[2], a[3]]]) * PERIODIC_TABLE[a[0]]["mass"])
        centro *= (1 / M)
        return centro

    def get_center_atom(self):
        """ Metodo para calcular el atomo central de la molecula. """
        distancias = []
        centro = self.get_center_of_mass()
        for i in xrange(len(self.atoms)):
            temp = self.atoms[i].coords - centro
            distancias.append([i, self.atoms[i].element, temp.norm()])
        distancias.sort(key=lambda s: s[2])
        return distancias[0]

    def get_center(self):
        """ Metodo para encontrar el centro geometrico de la molecula. """
        atms = len(self.atoms)
        centro = Matrix([[0,0,0]])
        for i in xrange(atms):
            centro = centro + self.atoms[i].coords
        centro *= (1.0/atms)
        return centro

    def get_limits(self):
        """ Metodo para obtener los extremos (en tamanio)
        de la molecula. """
        coords = self.get_coords()
        fin = len(coords) - 1
        coords.sort(key=lambda s: s[1])
        min_x, max_x = coords[0][1], coords[fin][1]
        coords.sort(key=lambda s: s[2])
        min_y, max_y = coords[0][2], coords[fin][2]
        coords.sort(key=lambda s: s[3])
        min_z, max_z = coords[0][3], coords[fin][3]
        return [min_x, max_x, min_y, max_y, min_z, max_z]

    def select(self, *selection):
        l = len(self.atoms)
        ctrl = []
        try:
            s = selection[0]
            if type(s) == int and s > 0 and s < l:
                ctrl.append(s)
            elif type(s) == str and s == "all":
                ctrl = range(l)
            elif type(s) == list:
                ctrl = s
            else:
                raise TypeError
            for i in ctrl:
                self.atoms[i].flag = not self.atoms[i].flag
            return True
        except:
            raise TypeError, "Molecule: The selection could not be made."

    def move(self, vector):
        try:
            vector = Matrix(vector)
        except:
            raise TypeError, "Molecule: The given value is not a vector."
        for a in self.atoms:
            temp = a.coords + vector
            a.set_coordinates(temp[0][0], temp[0][1], temp[0][2])
        return True

    def rotate(self, center, angle, vector):
        try:
            vector = Matrix(vector)
            center = Matrix(center)
        except:
            raise TypeError, "Molecule: The given value is not a vector."
        R = rotation_mat(angle, vector)
        self.move(-center)
        ctrl = 0
        for a in self.atoms:
            if a.flag:
                ctrl += 1
                a.coords = R * a.coords
        self.move(center)
        if ctrl == 0:
            raise RuntimeError, "Molecule: No atom was moved!"
        elif ctrl > len(self.atoms):
            raise IndexError, "Molecule: Too many atoms were moved."
        else:
            return True

    def make_grid(self, mesh=0.5, extent=4, charge=False, lims=False):
        """ Funcion para construir una red de puntos sensibles al
        campo electrico de la molecula. """
        if not lims:
            limits = self.get_limits()
            for i in xrange(6):
                if i % 2 == 0:
                    limits[i] -= extent
                else:
                    limits[i] += extent
        else:
            limits = [l for l in lims]
        if (type(mesh) == float) or (type(mesh) == int):
            m = [mesh, mesh, mesh]
        else:
            m = mesh
        d_x = int((limits[1] - limits[0]) / m[0])
        d_y = int((limits[3] - limits[2]) / m[1])
        d_z = int((limits[5] - limits[4]) / m[2])
        grid = [0] * d_x * d_y * d_z
        coords = self.get_coords()
        for i in xrange(d_x):
            for j in xrange(d_y):
                for k in xrange(d_z):
                    x = limits[0] + i * m[0]
                    y = limits[2] + j * m[1]
                    z = limits[4] + k * m[2]
                    if charge:
                        c = feel_field(x, y, z, coords)
                        grid[i*d_y*d_z + j*d_z + k] = (x, y, z, c)
                    else:
                        grid[i*d_y*d_z + j*d_z + k] = (x, y, z)
        return grid

    def make_bond_grid(self, bond, angle_mesh=math.pi/6, mesh=0.5, charge=False):
        """ Funcion para construir una red de puntos alrededor de un enlace,
        sensibles al campo electrico de la molecula. """
        pass

# Funciones especiales del paquete

# Metodo para calcular el campo electrico en un punto determinado de la molecula.
def feel_field(x, y, z, coords):
    vector = [0, 0, 0, 0]
    for c in coords:
        tvec = [0, 0, 0]
        tvec[0] = ((c[1] - x)**2)**0.5
        tvec[1] = ((c[2] - y)**2)**0.5
        tvec[2] = ((c[3] - z)**2)**0.5
        r = (tvec[0] ** 2 + tvec[1] ** 2 + tvec[2] ** 2) ** 0.5
        E = K_ELE * c[4] / r
        tvec = [tvec[i] * E / r for i in xrange(3)]
        vector[0] += tvec[0]
        vector[1] += tvec[1]
        vector[2] += tvec[2]
        vector[3] += E
    return vector

# Re-estructurador de enlaces
def rebonder(info):
    temp = []
    for i in info[1:]:
        try:
            temp.append(int(i))
        except Exception as e:
            if i == 'ar':
                temp.append(8)
            elif i == 'am':
                temp.append(9)
            else:
                raise TypeError, e
    temp[0] -= 1
    temp[1] -= 1
    return temp

# Modificador de elementos para busqueda
def elementor(a):
    if a in PERIODIC_TABLE.keys():
        return a
    else:
        if a[:2] in PERIODIC_TABLE.keys():
            return a[:2]
        else:
            if a[0] in PERIODIC_TABLE.keys():
                return a[0]
            else:
                msg = "The symbol {0} is not an element in the Periodic Table.".format(a)
                raise ValueError, msg

# Abrir molecula en formato mol2
def mol_from_mol2(smol):
    smol = str(smol).split('@')
    try:
        smol.remove('')
    except ValueError as e:
        pass
    ats = []
    bds = []
    for l in smol:
        if '<TRIPOS>ATOM' in l:
            lines = l.split('\n')[1:]
            lines = [line.split() for line in lines]
            ats = [Atom(element=elementor(line[1]), x=float(line[2]), y=float(line[3]),
                        z=float(line[4]), charge=float(line[-1])) for line in lines if line != []]
        elif '<TRIPOS>BOND' in l:
            lines = l.split('\n')[1:]
            bds = [rebonder(line.split()) for line in lines if len(line) != 0]
    new = Molecule()
    new.add_atoms(ats)
    new.set_bonds(bds)
    return new

def make_cube(molecule, mesh=0.5, fname='mol'):
    lims = molecule.get_limits()
    extent = 3.7
    for h in xrange(6):
        if h % 2 == 0:
            lims[h] -= extent
        else:
            lims[h] += extent
    D = [lims[1] - lims[0], lims[3] - lims[2], lims[5] - lims[4]]
    d = [int(r / mesh) for r in D]
    grid = molecule.make_grid(mesh=mesh, charge=True, lims=lims)
    cg = [list(g[:-1]) + [g[3][3]] for g in grid]
    cl = ['Cube file build with PyChemT by Rony J. Letona\n',
          'Grid built using electrostatics insted of MOs\n',
          '\t{0}\t{1:12.6f}\t{2:12.6f}\t{3:12.6f}\n'.format(len(molecule.atoms), lims[0]*BOHR,
                                                            lims[2]*BOHR, lims[4]*BOHR),
          '\t{0}\t{1:12.6f}\t{2:12.6f}\t{3:12.6f}\n'.format(d[0], BOHR*D[0]/(d[0]-1), 0.0, 0.0),
          '\t{0}\t{1:12.6f}\t{2:12.6f}\t{3:12.6f}\n'.format(d[1], 0.0, BOHR*D[1]/(d[1]-1), 0.0),
          '\t{0}\t{1:12.6f}\t{2:12.6f}\t{3:12.6f}\n'.format(d[2], 0.0, 0.0, BOHR*D[2]/(d[2]-1))]
    for i in xrange(len(molecule.atoms)):
        a = molecule.atoms[i]
        n = PERIODIC_TABLE[a.element]['number']
        crg = a.charge
        ax, ay, az = a.coords[0]
        cl.append('\t{0}\t{1:12.6f}\t{2:12.6f}\t{3:12.6f}\t{4:12.6f}\n'.format(n,
                                                crg, ax*BOHR, ay*BOHR, az*BOHR))
    temp = ''
    for j in range(d[0]):
        for k in range(d[1]):
            for l in range(d[2]):
                temp += '{0:13.5e}\t'.format(cg[(j*d[1] + k)*d[2] + l][3] * 10**9)
                if ((j*d[1] + k)*d[2] + l + 1) % 6 == 0:
                    temp += '\n'
            temp += '\n'
    cl.append(temp)
    with open(fname + '.cube', 'w') as f:
        f.writelines(cl)

# Matriz de rotacion para vectores
def rotation_mat(a=math.pi, u=Matrix([[1, 1, 1]])):
    R = [[math.cos(a) + u[0][0] ** 2 * (1 - math.cos(a)),
          u[0][0] * u[0][1] * (1 - math.cos(a)) + u[0][2] * math.sin(a),
          u[0][0] * u[0][2] * (1 - math.cos(a)) - u[0][1] * math.sin(a)],
         [u[0][0] * u[0][1] * (1 - math.cos(a)) - u[0][2] * math.sin(a),
          math.cos(a) + u[0][1] ** 2 * (1 - math.cos(a)),
          u[0][1] * u[0][2] * (1 - math.cos(a)) + u[0][0] * math.sin(a)],
         [u[0][0] * u[0][2] * (1 - math.cos(a)) + u[0][1] * math.sin(a),
          u[0][1] * u[0][2] * (1 - math.cos(a)) - u[0][0] * math.sin(a),
          math.cos(a) + u[0][2] ** 2 * (1 - math.cos(a))]]
    return Matrix(R)
