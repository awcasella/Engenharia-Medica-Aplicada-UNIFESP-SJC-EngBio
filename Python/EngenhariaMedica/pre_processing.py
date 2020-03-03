import numpy as np
from scipy.signal import welch
from scipy.stats import skew, kurtosis
from scipy.interpolate import Rbf
from itertools import permutations, combinations
import matplotlib.pyplot as plt 


def normalizacao(VETOR, metodo='std', r=1):
	""" Normalizes the values of a single feature.

	INPUT:
	- VETOR: valores de uma caracteristica calculada em N padroes e C classes,
	  Vetor com dimensao 1 x N*C
	- metodo ='std' : normalizacao linear (padrao)
	         = 'mmx': limitada entre -1 e 1
	         = 'sfm': rescala nao linear no intervalo 0 a 1
	- r = parametro do metodo sfm (padrao =1)
	
	OUTPUT:
	- VETORNORM = 1 x N*C: vetor com os valores normalizados da caracteristica 
	
	"""
	
	M, N = VETOR.shape
	VETOR = VETOR.reshape(1, M*N)
	if metodo == 'std':
		VETORNORM = VETOR-VETOR.mean()
		VETORNORM = VETORNORM/(VETOR.std())
	elif metodo == 'mmx':
		VETORNORM = 2*VETOR/(max(VETOR)-min(VETOR))
		VETORNORM = VETORNORM-(min(min(VETORNORM))+1)
	elif metodo == 'sfm':
		Y = VETOR-VETOR.mean()
		Y = Y/(r*VETOR.std())
		VETORNORM = 1/(1+np.exp(-Y))
	else:
		raise AttributeError("Unknown method, but don't get sad, not everything is made of roses.")
	VETORNORM = VETORNORM.reshape(M, N)

	return VETORNORM


def rmoutliers(data, p=3):
	""" Remove outliers de um conjunto de valores utilizando como limiar um numero p de desvios padroes em relacao a mediana.

	INPUT:
	- x = vector of data of a single feature.
	- p = number of standard deviations to be used.
	
	OUTPUT
	- data: data without outliers
	- outliers: outliers detected
	- indexes: indexes of outliers

	"""

	inferior = np.where(data<np.median(data) - p*data.std())[0]
	superior = np.where(data>np.median(data) + p*data.std())[0]
	indexes = np.union1d(inferior, superior)
	outliers = data[indexes]
	data = np.delete(data, indexes)

	return data, outliers, indexes


