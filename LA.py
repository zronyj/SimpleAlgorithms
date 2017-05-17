# Clase matriz para contar con un sistema en el que se pueda
# trabajar vectores.


class Matrix(list):
    def __repr__(self):
        """ Metodo para representar la matriz de una manera facil de entender. """
        c = len(self)
        f = len(self[0])
        out = "\t"
        for y in xrange(f):
            for x in xrange(c):
                out += str(round(self[x][y], 5)) + "\t"
            out += "\n\t"
        return out

    def __add__(self, other):
        """ Metodo para la suma de matrices. """
        if self.check("add", other):
            new = [map(lambda r, s: r + s, i, j) for i, j in zip(self, other)]
            return Matrix(new)
        else:
            raise TypeError, "Matrix: Dimensions do not match."

    def __sub__(self, other):
        """ Metodo para la resta de matrices. """
        if self.check("add", other):
            new = [map(lambda r, s: r - s, i, j) for i, j in zip(self, other)]
            return Matrix(new)
        else:
            raise TypeError, "Matrix: Dimensions do not match."

    def __neg__(self):
        """ Metodo para calcular la forma negativa de matrices. """
        new = [map(lambda r: -1 * r, i) for i in self]
        return Matrix(new)

    def __mul__(self, other):
        """ Metodo para la multiplicacion de matrices. """
        if type(other) == int or type(other) == float:
            new = [map(lambda r: other * r, i) for i in self]
            return Matrix(new)
        else:
            if self.check("mul", other):
                tra = self.transpose()
                fa = len(tra)
                cb = len(other)
                new = [[0] * fa for t in xrange(cb)]
                for f in xrange(fa):
                    for c in xrange(cb):
                        new[c][f] = sum(map(lambda r, s: r * s, tra[f], other[c]))
                return Matrix(new)
            else:
                raise IndexError, "Matrix: Dimensions are not a match to multiply."

    def __pow__(self, other):
        """ Metodo para operar potencias de matrices. """
        if type(other) != int:
            raise TypeError, "Matrix: Can not operate fractional power on matrices."
        else:
            if other >= 1:
                new = Matrix(self * 1)
                for i in xrange(other - 1):
                    new = new * self
                return new
            elif other == 0:
                a = len(self)
                b = len(self[0])
                new = [[0] * b for t in xrange(a)]
                if a >= b:
                    c = b
                else:
                    c = a
                for i in xrange(c):
                    new[i][i] = 1
                return Matrix(new)
            else:
                raise ValueError, "Matrix: Can not operate negative power on matrices."

    def check(self, op, other):
        """ Metodo para revisar si dos matrices son operables. """
        if op == "add" or op == "sub":
            ia, ja = self.dim()
            ib, jb = other.dim()
            if ia == ib and ja == jb:
                return True
            else:
                return False
        elif op == "mul":
            fila = len(self)
            columna = len(other[0])
            if fila == columna:
                return True
            else:
                return False
        else:
            raise ValueError, "Matrix: There is no such operation."

    def diag(self):
        """ Metodo para extraer la diagonal principal de una matriz. """
        if len(self) > len(self[0]):
            d = len(self)
        else:
            d = len(self[0])
        temp = [0] * d
        for i in xrange(d):
            temp[i] = self[i][i]
        return temp

    def row(self, i):
        """ Metodo para extraer una fila determinada de la matriz. """
        return self.transpose()[i]

    def col(self, i):
        """ Metodo para extraer una columna determinada de la matriz. """
        return self[i]

    def round_m(self, ndigits=10):
        temp = Matrix([[round(j, ndigits) for j in i] for i in self])
        return temp

    def transpose(self):
        """ Metodo para transponer la matriz. """
        temp = zip(*list(self))
        new = [list(e) for e in temp]
        return Matrix(new)

    def t(self):
        """ Metodo rapido para transponer la matriz. """
        return self.transpose()

    def trace(self):
        """ Metodo para calcular la traza de una matriz. """
        return sum(self.diag())

    def submat(self, i, j):
        """ Metodo para generar una submatriz (una menor) de la matriz original. """
        c = len(self)
        r = len(self[0])
        new = Matrix([[0] * (r - 1) for t in xrange(c - 1)])
        for k in xrange(i):
            for l in xrange(j):
                new[k][l] = self[k][l]
            for m in xrange(j + 1, r):
                new[k][m - 1] = self[k][m]
        for n in xrange(i + 1, c):
            for l in xrange(j):
                new[n - 1][l] = self[n][l]
            for m in xrange(j + 1, r):
                new[n - 1][m - 1] = self[n][m]
        return new

    def cofactor_matrix(self):
        """ Metodo para calcular la matriz de cofactores. """
        c = len(self)
        r = len(self[0])
        new = Matrix([[0] * r for t in xrange(c)])
        for x in xrange(c):
            for y in xrange(r):
                new[x][y] = (-1) ** (x + y) * self.submat(x, y).det()
        return new

    def det(self):
        """ Metodo para calcular la determinante de la matriz. """
        r = len(self[0])
        c = len(self)
        if r == c:
            if r == 2:
                return self[0][0] * self[1][1] - self[1][0] * self[0][1]
            elif r > 2:
                d = 0
                for i in xrange(c):
                    d += (-1) ** i * self[i][0] * self.submat(i, 0).det()
                return d
            else:
                raise IndexError, "Matrix: Error in dimensions."
        else:
            raise IndexError, "Matrix: Dimensions do not match."

    def dim(self):
        """ Metodo para indicar las dimensiones de la matriz: fila, columna. """
        return len(self[0]), len(self)

    def adj(self):
        """ Metodo para calcular la adjunta de una matriz. """
        return Matrix(self.cofactor_matrix().transpose())

    def inv(self):
        """ Metodo para calcular la inversa de una matriz. """
        r = len(self[0])
        c = len(self)
        if r == c:
            if r == 2:
                return Matrix([[self[1][1], -1 * self[0][1]],
                               [-1 * self[1][0], self[0][0]]]) * (1.0 / self.det())
            elif r > 2:
                return Matrix(self.adj() * (1.0 / self.det()))
            else:
                raise IndexError, "Matrix: Error in dimensions."
        else:
            raise IndexError, "Matrix: Dimensions do not match."

    def pinv(self):
        """ Metodo para calcular la inversa de una matriz.
        utilizando el metodo de Moore-Penrose y descomposicion QR. """
        Q, R = self.QR()
        return (R.transpose() * R).inv() * self.transpose()

    def v_proj(self, other):
        """ Metodo para proyectar un vector sobre otro. """
        if len(self) == 1 and len(other) == 1 and self.check("add", other):
            up = self.t() * other
            down = other.t() * other
            coef = float(up[0][0]) / down[0][0]
            return other * coef
        else:
            raise TypeError, "Matrix: Can only operate on vectors."

    def norm(self):
        """ Metodo para calcular la norma de una matriz. """
        temp = sum([sum(map(lambda r: r ** 2, i)) for i in self])
        return temp ** 0.5

    def unitary(self):
        """ Metodo para calcular un vector unitario. """
        n = self.norm()
        if n == 0: n = 1
        return self * (1.0 / self.norm())

    def QR(self):
        """ Metodo para calcular una descomposicion QR. """
        x = len(self)
        Q = []
        for i in xrange(x):
            v = self.col(i)
            if len(Q) == 0:
                ein = Matrix([v]).unitary()[0]
                Q.append(ein)
            else:
                u = Matrix([v])
                for j in Q:
                    w = Matrix([j])
                    u -= u.v_proj(w)
                ein = u.unitary()[0]
                Q.append(ein)
        Q = Matrix(Q)
        R = (Q.t() * self)
        return Q, R

    def Eigen(self):
        """ Metodo para calcular valores y vectores propios mediante
        descomposicion QR. """
        M = []
        M.append(self)
        Q, R = M[0].QR()
        M.append(R * Q)
        EV = Q
        while Matrix([(M[-1] - M[-2]).diag()]).norm() > 1e-10:
            Q, R = M[-1].QR()
            M.append(R * Q)
            EV *= Q
        return M[-1].round_m(6), EV

    def SVD(self):
        """ Metodo para hallar la descomposicion de valores singulares. """
        ctc = self.transpose() * self
        S, V = ctc.Eigen()
        eigen_vals = S.diag()
        x = len(S)
        S *= 0
        for s in xrange(x):
            S[s][s] = eigen_vals[s] ** 0.5
        U = self * V * S.inv()
        return U, S, V