import numpy as np
import scipy.linalg as sla


class Interpolator(object):

	"""Cubic Hermite Spline Interpolation"""

	def __init__(self, x, y):

		"""
		Public variables: locations (x) and values (y) of knots, 
		derivatives (d) at knots, and quadratic (c) and cubic (b) 
		term coefficients
		"""

		self.x = x
		self.y = y
		self.__h = x[1:] - x[:-1]
		self.__delta = (y[1:] - y[:-1]) / self.__h
		self.d = self.__set_derivatives(self.__h, self.__delta)
		self.c = (3.0 * self.__delta - 2.0 * self.d[0: -1] - \
				  self.d[1:]) / self.__h
		self.b = (self.d[0: -1] - 2.0 * self.__delta + \
				  self.d[1:]) / self.__h ** 2

	def __set_derivatives(self, h, delta):

		"""
		Calculate derivatives at each knot, using not-a-knot 
		formalism (see www.mathworks.com/content/dam/mathworks/
		mathworks-dot-com/moler/interp.pdf)
		"""

		# dimension compnents
		n_diffs = len(h) + 1
		r_vec = np.zeros(n_diffs)
		a_mat = np.zeros((n_diffs, n_diffs))

		# form delta vector
		r_vec[0] = ((h[0] + 2.0 * (h[0] + h[1])) * h[1] * delta[0] + \
					h[0] ** 2 * delta[1]) / (h[0] + h[1])
		r_vec[1: -1] = 3.0 * (h[1:] * delta[0: -1] + \
							  h[0: -1] * delta[1:])
		r_vec[-1] = (h[-1] ** 2 * delta[-2] + \
					 (2.0 * (h[-2] + h[-1]) + h[-1]) * \
					 h[-2] * delta[-1]) / (h[-2] + h[-1])
		
		'''
		print h
		print delta
		print r_vec
		print 3 * (delta[0] * 5 + delta[1]) / 6
		print 3 * (delta[0:-2] + delta[1:-1])
		print 3 * (delta[-1] * 5 + delta[-2]) / 6
		'''

		# form dense A matrix
		a_mat[0, 0] = h[1]
		a_mat[0, 1] = h[1] + h[0]
		for i in range(1, n_diffs - 1):
			a_mat[i, i-1] = h[i]
			a_mat[i, i] = 2.0 * (h[i-1] + h[i])
			a_mat[i, i+1] = h[i-1]
		a_mat[-1, -2] = h[-1] + h[-2]
		a_mat[-1, -1] = h[-2]
		
		#print a_mat

		# solve system and return
		return np.linalg.solve(a_mat, r_vec)

	def __set_derivatives_banded(self, h, delta):

		"""
		Calculate derivatives at each knot, using not-a-knot 
		formalism and forming a banded (not dense) matrix
		"""

		# dimension components
		n_diffs = len(h) + 1
		r_vec = np.zeros(n_diffs)
		a_mat = np.zeros((n_diffs, n_diffs))

		# form delta vector
		r_vec[0] = ((h[0] + 2.0 * (h[0] + h[1])) * h[1] * delta[0] + \
					h[0] ** 2 * delta[1]) / (h[0] + h[1])
		r_vec[1: -1] = 3.0 * (h[1:] * delta[0: -1] + \
							  h[0: -1] * delta[1:])
		r_vec[-1] = (h[-1] ** 2 * delta[-2] + \
					 (2.0 * (h[-2] + h[-1]) + h[-1]) * \
					 h[-2] * delta[-1]) / (h[-2] + h[-1])

		# form banded A matrix
		a_mat_ud = np.concatenate(([0], [h[1] + h[0]], h[0: -1]))
		a_mat_d = np.concatenate(([h[1]], 2.0 * (h[0: -1] + h[1:]), [h[-2]]))
		a_mat_ld = np.concatenate((h[1:], [h[-1] + h[-2]], [0]))
		a_mat_band = np.matrix([a_mat_ud, a_mat_d, a_mat_ld])

		# solve system and return
		return sla.solve_banded((1, 1), a_mat_band, r_vec, \
								overwrite_ab=True, \
								overwrite_b=True)

	def __interpolate_scalar(self, x_int, der=False):

		"""Interpolate an individual point"""

		ind = [i for i, check in enumerate(x_int < self.x) if check][0] - 1
		s = x_int - self.x[ind]
		if der:
			y_int = self.d[ind] + \
					s * (2.0 * self.c[ind] + \
						 s * 3.0 * self.b[ind])
		else:
			y_int = self.y[ind] + \
					s * (self.d[ind] + \
						 s * (self.c[ind] + \
							  s * self.b[ind]))
		return y_int

	def interpolate(self, x_int):

		"""Public interpolation method, accepts scalars or 1D arrays"""

		if isinstance(x_int, np.ndarray):
			if (x_int < self.x[0]).any() or \
			   (x_int > self.x[-1]).any():
				raise ValueError('Requested interpolation out of bounds!')
			n = len(x_int)
			y_int = np.zeros(n)
			for i in range(n):
				y_int[i] = self.__interpolate_scalar(x_int[i])
		else:
			if x_int < self.x[0] or x_int > self.x[-1]:
				raise ValueError('Requested interpolation out of bounds!')
			y_int = self.__interpolate_scalar(x_int)
		return y_int
