import numpy as np
from scipy.signal import welch
from scipy.stats import skew, kurtosis
from scipy.interpolate import Rbf
from itertools import permutations, combinations
import matplotlib.pyplot as plt 


def bayesian(medias, covariancias, priors, x):
	"""Classifier based on gaussian distribution of the features
	INPUTS
	- medias: 2D matrix with each column containing the mean of each class.
	- covariancias: 3D Matrix of 2D covariance matrixes of each class.
	- priors: wights of each class.
	- x: New pattern to be classified.

	OUTPUTS
	- P: Probability of new pattern belong to each class.
	- onde: Class in which new pattern is more likely to belong.
	"""

	M = covariancias.shape[0] # Number of classes
	
	lh = np.zeros(M)
	for n in range(0, M):
		Mcov = covariancias[n, :, :]
		mu = medias[:, n]
		norm = 1/np.sqrt(((2*np.pi)**L)*np.linalg.det(Mcov))
		lh[n] = norm*np.exp(-0.5*(x - mu).dot(np.linalg.inv(Mcov)).dot((x - mu).T))

	mlh = lh.dot(priors)
	P = lh*priors/mlh
	onde = int(np.where(P ==  max(P)))
	P = tuple(P)
	return P, onde


def minimumDistance(medias, Mcov, x):
	"""Minunimum distance classifier based on eucliidean and Mahallanobis distance
	INPUTS
	- medias: 2D matrix with each column containing the mean of each class.
	- Mcov: 2D Covariance matrix of whole space
	- x: New pattern to be classified.

	OUTPUTS
	- mahalanobis: Mahalanobis distance of new pattern to each class.
	- mahalanobisClass: Classification of new patter accordingly to mahalanobis distance.
	- euclideana: Euclidean distance of new pattern to each class.
	- euclideanaClass: Classification of new patter accordingly to euclidean distance.
	"""
	M = medias.shape[1]
	euclideana = np.zeros(M)
	mahalanobis = np.zeros(M)
	for n in range(0, M):
		mu = medias[:, n]
		euclideana[n] = (x - mu).dot((x - mu).T)
		mahalanobis[n] = (x - mu).dot(np.linalg.inv(Mcov)).dot((x - mu).T)
	euclideanaClass = int(np.where(euclideana == min(euclideana)))
	mahalanobisClass = int(np.where(mahalanobis == min(mahalanobis)))
	
	mahalanobis = tuple(mahalanobis)
	euclideana = tuple(euclideana)

	return mahalanobis, mahalanobisClass, euclideana, euclideanaClass


def perceptron(classe1, classe2, rho=0.1, niter=1000):
	"""Binary Perceptron classifier 
	INPUTS
	- classe1: 2D numpy matrix of class1
	- classe2: 2D numpy matrix of class1
	- rho: Weight of error, indicates how much is the step during each iteration.
	- niter: Number of maximum iterations to perform.
	
	OUTPUTS
	- w: Perceptron hyperplane.
	- iter: Number of iterations it took to calculate w
	"""

	J = lambda p, w, dados, Y: sum(p[Y]*w.dot(dados[:, Y]))
	L, N1 = classe1.shape
	L2, N2 = classe2.shape
	iter = niter

	c = np.ones([1, N1+N2])
	c[0, N1:] = -1

	if L != L2:
		print('Precisa que as classes tenham o mesmo numero de caracteristicas.')
		return
	dados = np.ones([L+1, N1+N2])
	dados[:-1, :N1] = classe1
	dados[:-1, N1:] = classe2

	w = np.random.randn(L+1)
	p = np.sign(w.dot(dados))

	inutil, Y = np.where(c != p)
	plt.figure()
	e = []
	erroAntes = np.zeros(L+1)
	erroSum = np.zeros(L+1)
	Kp = 0.05
	Ki = 0
	Kd = 0
	for n in range(0, niter):
		erro = sum((p[Y]*dados[:, Y]).T)
		P = rho*erro
		I = Ki*(erroSum + erro)
		D = Kd*(erro - erroAntes)/0.4
		PD = P + I + D
		w = w - PD
		erroAntes = erro
		p = np.sign(w.dot(dados))
		e.append(J(p, w, dados, Y))
		inutil, Y = np.where(c != p)
		if np.where(c!=p)[0].shape[0] == 0:
			iter = n+1
			break
	plt.plot(e)
	plt.title('Erro perceptron')
	plt.show()
	return w, iter


def pocket(classe1, classe2, rho=0.1, niter=1000):
	"""Binary Perceptron Pocket classifier
	INPUTS
	- classe1: 2D numpy matrix of class1
	- classe2: 2D numpy matrix of class1
	- rho: Weight of error, indicates how much is the step during each iteration.
	- niter: Number of maximum iterations to perform.
	
	OUTPUTS
	- w0: Perceptron hyperplane.
	- iter: Number of iterations it took to calculate w
	"""

	J = lambda p, w, dados, Y: sum(p[Y]*w.dot(dados[:, Y]))
	L, N1 = classe1.shape
	L2, N2 = classe2.shape
	iter = niter

	c = np.ones([1, N1+N2])
	c[0, N1:] = -1

	if L != L2:
		print('Precisa que as classes tenham o mesmo numero de caracteristicas.')
		return
	dados = np.ones([L+1, N1+N2])
	dados[:-1, :N1] = classe1
	dados[:-1, N1:] = classe2

	w0 = np.random.randn(L+1)
	p = np.sign(w0.dot(dados))

	inutil, Y = np.where(c != p)
	e = []
	erroAntes = np.zeros(L+1)
	erroSum = np.zeros(L+1)
	Kp = rho
	Ki = 0
	Kd = 0
	for n in range(0, niter):
		erro = sum((p[Y]*dados[:, Y]).T)
		P = rho*erro
		I = Ki*(erroSum + erro)
		D = Kd*(erro - erroAntes)/0.4
		PD = P + I + D
		wp = w0 - PD
		erroAntes = erro
		e.append(J(p, wp, dados, Y))
		if J(p, wp, dados, Y) < J(p, w0, dados, Y):
			w0 = wp
		
		p = np.sign(w0.dot(dados))
		inutil, Y = np.where(c != p)
		if np.where(c!=p)[0].shape[0] == 0:
			iter = n+1
			break
	plt.figure()
	plt.plot(e)
	plt.title('Erro pocket')
	plt.show()
	return w0, iter
	

def mle(dados, classes):
	M = int(classes[0][-1] + 1) # Numero de classes
	L = int(dados.shape[0]) # Numero de caracteristicas
	media = np.zeros([L, M])
	covariancias = np.zeros([M, L, L])
	for n in range(0, M): # Estima a media e covariancia de cada classe
		media[:, n] = np.average(dados,axis=1)
		covariancias[n, :, :] = np.cov(dados)

	return media, covariancias


def ls(x1, x2, alpha):
	y = np.concatenate((x1, x2), axis=1)
	uns = np.ones([1, y.shape[1]])
	y = np.concatenate((y, uns), axis=0)
	L,N = y.shape
	uns = np.ones([1, x1.shape[1]])
	menosUns = -1*np.ones([1, x2.shape[1]])
	c = np.concatenate((uns, menosUns), axis=1)
	termo1 = np.zeros([L,L]);
	termo2 = np.zeros([L,1]);
	for k in range(0, N):
		termo1 += y[:, k].dot(y[:, k].T) + alpha*np.eye(L)
		termo2 += c[0, k]*np.reshape(y[:, k], (-1, 1))
	w = np.linalg.inv(termo1).dot(termo2).T

	return w

def scatter(classes, K, axis=0):
	C = len(classes)
	
	if axis == 1:
		for n in range(0, C):
			classes[n] = classes[n].T
	
	L = classes[0].shape[0]
	Nc = np.zeros([C])
	dados = classes[0]
	for n in range(0, C):
		c = classes[n]
		Nc[n] = c.shape[1]
		if n > 0:
			dados = np.concatenate((dados, c), axis=1)

	N = sum(Nc)
	P = Nc/N
	
	Sw = np.zeros([K, K])
	Sm = np.zeros([K, K])
	Sb = np.zeros([K, K])

	for subset in combinations(range(0, L), K):
		for n in range(0, C):
			c = classes[n]
			Mcov = np.cov(c[subset, :], ddof=0)
			Sw += P[n]*Mcov

		Sm += np.cov(dados[subset, :], ddof=0)

	Sb = Sm - Sw
	return Sw, Sb

def fda(classes, n):
	"""FDA 
	
	INPUTS:
	- classes: List of 2D numpy matrixes of classes.
	- n: Number of features in the new space.
	
	OUTPUTS:
	- Y: New optmized space.
	"""

	L = classes[0].shape[1] # Number of features
	Sw, Sb = scatter(classes, L)
	eigval, eigvec = np.linalg.eig(np.linalg.inv(Sw).dot(Sb))

	eigvec, eigval = autoordenar(eigvec, eigval)

	C = len(classes)                               # Number of classes.
	A = eigvec[:, 0:C-1].T

	X = classes[0]

	for n in range(1, C):
		X = np.concatenate((X, classes[n]), axis=1)


	Y = A.dot(X)

	return Y
