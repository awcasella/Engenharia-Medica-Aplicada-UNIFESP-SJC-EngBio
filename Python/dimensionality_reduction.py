import numpy as np
from scipy.signal import welch
from scipy.stats import skew, kurtosis
from scipy.interpolate import Rbf
from itertools import permutations, combinations
import matplotlib.pyplot as plt 



def aula25_GerandoDadosGaussianos(medias, covariancias, N, priors=None, plotar=1, seed=0):
	""" Generates gaussian N-point data based on the mean and covariances, and priors
	
	INPUTS
	- medias: Numpy array of mean of each class in its rows.
	- covariancias: List of each class' covariance matrix.
	- N: Number of patterns to be generated.
	- priors: Weights for each class.
	- plotar: Flag to plot the data or not.
	- seed: Seed for generating pseudorandom numbers.

	OUTPUTS
	- dadossim: Simulated data, each row containing a pattern and each column containing a feature.
	- classessim: Numpy array of values of class of each pattern generated.
	"""
	
	L = medias.shape[0]
	M = medias.shape[1]

	if priors is None:
		priors = np.ones([medias.shape[1], 1])/medias.shape[1]
	
	Ni = np.zeros([M+1, 1])
	for i in range(1, M+1):
		n = i-1
		Ni[i] = np.around(priors[n]*N)

	np.random.seed(seed)
	dadossim = np.zeros([int(sum(Ni)), L])
	classessim = np.zeros([1, int(sum(Ni))])

	for i in range(1, M+1):
		n = i-1
		dadossim[int(sum(Ni[0:i])):int(sum(Ni[0:i+1])), :] = np.random.multivariate_normal(medias[:,n], covariancias[n,:,:], int(Ni[i]))
		classessim[:, int(sum(Ni[0:i])):int(sum(Ni[0:i+1]))] = n

	return dadossim, classessim


def autoordenar(eigvec,eigval):
	#eigvalord = sorted(eigval, reverse=True)
	idx = np.argsort(eigval)
	eigvec = eigvec.T
	eigvecord = eigvec[idx]
	eigvalord = eigval[idx]
	eigvecord = eigvecord[-1::-1]
	eigvalord = eigvalord[-1::-1]
	return eigvecord.T, eigvalord


def klt(X, n):
	"""Calculates the n principal components of X
	INPUTS
	- X: Data matrix, features are on the columns
	- n: Number of principal components to be calculated.

	OUTPUTS
	- eigvec: Eigenvectors which are a base for the new space generated.
	- eigval: Eigenvalues related to the eigenvectors.
	- Y: Data projected in the eigenvectors.
	- erro: error due to loss of least important principal components.
	"""
	X = (X - X.mean(axis=0)) # Remover a media de cada caracteristica para centralizar os dados
	N = X.shape[0]
	Mcov = X.T.dot(X)/N # Matriz de covariancia
	eigval, eigvec = np.linalg.eig(Mcov)
	eigvec, eigval = autoordenar(eigvec,eigval)
	A = eigvec.T # Autovetores nas linhas
	Y = A.dot(X.T).T
	erro = sum(eigval[n:])/sum(eigval)
	return eigvec, eigval, Y[:, :n], erro


