import math, random

class Particle(object):

	def __init__(self, lims, func):
		try:
			self.coords = []
			self.limits = lims
			self.eval = func
			dims = len(lims)
			for i in range(dims):
				lower, upper = lims[i][0], lims[i][1]
				diff = upper - lower
				if diff < 0:
					raise Exception("Limits are wrong.")
				self.coords.append(lower + diff * random.random())
			self.value = self.eval(self.coords)
		except Exception as err:
			print err
			return None

	def __repr__(self):
		output = "Value: " + str(self.value) + "\n"
		for i in range(len(self.coords)):
			output += "q(" + str(i+1) + ") = " + str(self.coords[i]) + " ... " + str(self.limits[i]) + "\n"
		return output

	def new_coordinates(self):
		ctrl = []
		for i in range(len(self.coords)):
			ctrl.append(-1 + 2 * random.random())
		return ctrl

	def update_coordinates(self, new_coords):
		if len(new_coords) == len(self.coords):
			for i in range(len(self.coords)):
				self.coords[i] = new_coords[i]
			self.value = self.eval(self.coords)
		else:
			raise IndexError("Dimensions of vectors are not the same.")

def vect_add(vo, vf):
	nv = []
	for i in range(len(vf)):
		nv.append(vf[i] + vo[i])
	return nv

def vect_sub(vo, vf):
	nv = []
	for i in range(len(vf)):
		nv.append(vf[i] - vo[i])
	return nv

def scalar_mul(s, v):
	nv = []
	for i in range(len(v)):
		nv.append(s * v[i])
	return nv

def simulated_annealing(parts, lim, tempi, steps, optfunc, temp_prog=0):
	particles = []
	temp_funt = []
	const1 = steps / math.log(tempi + 1.0)
	const2 = float(tempi) / steps ** 2
	const3 = - tempi / float(steps)
	temp_funt.append( lambda t: math.e ** (- (t - steps) / const1) - 1 )
	temp_funt.append( lambda t: const2 * (t - steps) ** 2 )
	temp_funt.append( lambda t: const3 * t + tempi )
	f = open("path.pth", "w")
	f.writelines("Path of the particles.\n")
	f.close()
	try:
		cooling = temp_funt[temp_prog]
	except IndexError as err:
		print err
		return None
	for i in range(parts):
		particles.append(Particle(lim, optfunc))
	for j in range(steps):
		f = open("path.pth", "a")
		for k in range(parts):
			f.writelines("p(" + str(k) + ") : " + str(particles[k].coords) + " - " + str(particles[k].value) + "\n")
			v = particles[k].new_coordinates()
			u = vect_add(particles[k].coords, scalar_mul(cooling(j) / tempi, v))
			control = 0
			for l in range(len(u)):
				if u[l] < lim[l][0] or u[l] > lim[l][1]:
					control += 1
			if control == 0:
				energ = optfunc(u)
				if energ < particles[k].value:
					particles[k].update_coordinates(u)
				elif energ >= particles[k].value:
					cte = math.e ** (- (energ - particles[k].value) * tempi / cooling(j))
					if cte >= random.random() and cte < 1:
						particles[k].update_coordinates(u)
		f.close()
	for m in particles:
		print "---------------------------\n", m
	return particles

fv = lambda v: (v[0] - 11.2)**2 + (v[1] + 15.3)**2

simulated_annealing(4, [[0,15],[-20,0]], 100, 100000, fv, 0)
