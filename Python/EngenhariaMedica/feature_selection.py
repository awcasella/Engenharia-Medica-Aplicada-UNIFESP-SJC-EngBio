import numpy as np
from scipy.signal import welch
from scipy.stats import skew, kurtosis
from scipy.interpolate import Rbf
from itertools import permutations, combinations
import matplotlib.pyplot as plt 


def rocMeBabe(classe1, classe2):
	"""Computes ROC curve and AUC
	
	INPUT
	- classe1: Numpy array with first class of a feature.
	- classe2: Numpy array with second class of a feature.

	OUTPUT
	- AUC: Area under the ROC curve.
	- ACC: Accuracy.
	- VP: True positive.
	- VN: True negative.
	- FP: False positive.
	- FN: False negative.
	"""

	if classe1[-1] > classe2[-1]:
		aux = classe1
		classe1 = classe2
		classe2 = aux

	dist = np.union1d(classe1, classe2)
	VP = np.zeros([len(dist)+1])
	VN = np.zeros([len(dist)+1])
	FP = np.zeros([len(dist)+1])
	FN = np.zeros([len(dist)+1])
	for n in range(1, len(dist)+1):
		if dist[n-1] in classe1 and dist[n-1] not in classe2:
			VP[n] = VP[n-1] + 1
			VN[n] = FN[n-1] + 1
			FP[n] = FP[n-1]
			FN[n] = FN[n-1]
		if dist[n-1] in classe1 and dist[n-1] in classe2:
			VP[n] = VP[n-1] + 1
			VN[n] = VN[n-1]
			FP[n] = FP[n-1]
			FN[n] = FN[n-1] + 1
		if dist[n-1] not in classe1 and dist[n-1] in classe2:
			VP[n] = VP[n-1]
			VN[n] = VN[n-1]
			FP[n] = FP[n-1] + 1
			FN[n] = FN[n-1] + 1
			
	VP, VN, FP, FN = VP/max(VP), VN/max(VN), FP/max(FP), FN/max(FN)
	
	AUC = sum(((VP[1:]+VP[:-1])/2)*np.diff(FP, n=1))
	ACC = sum(VP+VN)/sum(VP+VN+FP+FN)

	return AUC, ACC, VP, VN, FP, FN


def fdr(classe1, classe2):
	""" Computes FDR criteria for two classes, this tells how apart each class is of each other

	INPUT
	- classe1: Numpy array with first class of a feature.
	- classe2: Numpy array with second class of a feature.

	OUTPUT
	- FDR: Value calculated for fdr of both classes.
	"""
	
	#return ((classe1.mean()-classe2.mean())**2)/(classe1.var(ddof=1) + classe2.var(ddof=1))
	if len(classe1.shape) == 1:
		classe1 = np.array([classe1])
	if len(classe2.shape) == 1:
		classe2 = np.array([classe2])
	
	num = ((np.average(classe1, axis=1) - np.average(classe2, axis=1))**2)
	den = (np.var(classe1, axis=1, ddof=1) + np.var(classe2, axis=1, ddof=1))
	FDR = tuple(num/den)
	return FDR


def selecaoEscalar(Mcorr, criterios, N=0, a1=0.5, a2=0.5):
	""" Performs a scalar feature selection which orders all features individually,
	from the best to the worst to separate the classes.

	INPUTS
	- Mcorr: Correlation matrix of all features.
	- criterios: 
	- N: Number of best features to be returned.
	- a1: Weigth for criterios.
	- a2: Weight for Mcorr.

	OUTPUTS
	- ordem: Tuple with the order of features.
	- M: Tuple with criteria for each feature.
	"""
	L = Mcorr.shape[0]
	if len(criterios.shape) != 1:
		criterios = criterios[0]
	if N==0 or N > len(criterios):
		N = len(criterios)
		print('You either did not specify or you gave a number grater than the number of characteristics.')
		print('Function will return all {} characteristics.'.format(N))
	
	Mcorr = abs(Mcorr)
	ordem = []
	M = []

	ordem.append(int(np.where(criterios == max(criterios))[0]))
	M.append(criterios[int(ordem[0])])
	Mcorr[:, int(ordem[0])] = 1

	fator = np.zeros(N) 
	for n in range(1, N):
		index = np.linspace(0, L-1, L)
		fator = np.sum(Mcorr[tuple(ordem), :], axis=0)
		MK = a1*criterios - a2*fator/n
		MK = np.delete(MK, ordem)
		index = np.delete(index, ordem)
		M.append(max(MK))
		ordem.append(int(index[int(np.where(MK == max(MK))[0])]))

	ordem = tuple(ordem)
	M = tuple(M)
	return ordem, M


def exaustivosel(classes, K, criterio):
	""" Selects the best set of k features which separate the classes.

	INPUTS
	- classes: Python List of numpy matrixes of each class, row are patterns and columns are features
	- K: Number of features to be selected.
	- metodo: Type of method to be used: 'exaustivo', 'forward' or 'floating'
	- criterio: Criteria to be used to calculate the best set of features.

	OUTPUTS
	- ordem: Set of feature which were selected.
	- maxcriterio: Value of criteria for the order calculated.
	"""

	L = classes[0].shape[1] # numero de caracteristicas 
	M = len(classes) # numero de classes
	
	Nc = np.zeros(M) # Numero de padroes em cada classes
	dados = classes[0]
	for n in range(0, M):
		c = classes[n]
		Nc[n] = c.shape[0]
		if n > 0:
			dados = np.concatenate((dados, c), axis=0)

	N = sum(Nc) # Total de padroes
	Pc = Nc/N  # Prob de padroes em cada classe
	
	maxcriterio = -np.inf

	for subset in combinations(range(0, L), K):
		Sw = np.zeros([K, K])
		Sb = np.zeros([K, K])
		Sm = np.zeros([K, K])

		for n in range(0, M):
			c = classes[n]
			matriz = np.cov(c[:, subset].T, ddof=0)
			Sw += Pc[n]*matriz

		Sm = np.cov(dados[:, subset].T, ddof=0)
		Sb = Sm - Sw

		if criterio.upper() == 'J1':
			J1 = Sm.trace()/Sw.trace()
			if J1 > maxcriterio:
				maxcriterio = J1
				ordem = subset[:]
		elif criterio.upper() == 'J2':
			J2 = np.linalg.det(np.linalg.inv(Sw).dot(Sm))
			if J2 > maxcriterio:
				maxcriterio = J2
				ordem = subset[:]
		elif criterio.upper() == 'J3':
			J3 = (np.linalg.inv(Sw).dot(Sm)).trace()/K
			if J3 > maxcriterio:
				maxcriterio = J3
				ordem = subset[:]
	
	ordem = tuple(ordem)
	return ordem, maxcriterio


def forwardsel(classes, K, criterio):
	""" Selects the best set of k features which separate the classes.

	INPUTS
	- classes: Python List of numpy matrixes of each class, row are patterns and columns are features 
	- K: Number of features to be selected.
	- metodo: Type of method to be used: 'exaustivo', 'forward' or 'floating'
	- criterio: Criteria to be used to calculate the best set of features.

	OUTPUTS
	- ordem: Set of feature which were selected.
	- maxcriterio: Value of criteria for the order calculated.
	"""

	L = classes[0].shape[1]
	M = len(classes)

	Nc = np.zeros(M)
	dados = classes[0]
	for n in range(0, M):
		c = classes[n]
		Nc[n] = c.shape[0]
		if n > 0:
			dados = np.concatenate((dados, c), axis=0)

	Pc = Nc/sum(Nc)

	ordem = []
	for S in range(1, K+1):
		maxcriterio = -np.inf
		for subset in combinations(range(0, L), S):
			# Verificando se "subset" contem as caracteristicas salvas em ordem
			if tuple(np.intersect1d(subset,ordem)) == tuple(sorted(ordem)):
				nova = int(np.setdiff1d(subset,ordem)) # caso contenha, acho a nova que sera testada
			else:
				continue # caso nao contenha, pulo essa iteracao

			Sw = np.zeros([S, S])
			Sm = np.zeros([S, S])
			Sb = np.zeros([S, S])

			for n in range(0, M):
				c = classes[n]
				Sw += Pc[n]*np.cov(c[:, subset].T, ddof=0)

			Sm = np.cov(dados[:, subset].T, ddof=0)

			if criterio.upper() == 'J1':
				if S == 1:
					J1 = Sm/Sw
				else:
					J1 = Sm.trace()/Sw.trace()
				if J1 > maxcriterio:
					maxcriterio = J1
					vaiessa = nova
			elif criterio.upper() == 'J2':
				if S == 1:
					J2 = Sm/Sw
				else:
					J2 = np.linalg.det(np.linalg.inv(Sw).dot(Sm))
				if J2 > maxcriterio:
					maxcriterio = J2
					vaiessa = nova
			elif criterio.upper() == 'J3':
				if S == 1:
					J3 = Sm/Sw
				else:
					J3 = np.trace(np.linalg.inv(Sw).dot(Sm))/S
				if J3 > maxcriterio:
					maxcriterio = J3
					vaiessa = nova

		ordem.append(vaiessa)
	
	ordem = tuple(ordem)
	return ordem, maxcriterio


def floatingsel(classes, K, criterio):
	""" Selects the best set of k features which separate the classes.

	INPUTS
	- classes: Python List of numpy matrixes of each class, row are patterns and columns are features
	- K: Number of features to be selected.
	- metodo: Type of method to be used: 'exaustivo', 'forward' or 'floating'
	- criterio: Criteria to be used to calculate the best set of features.

	OUTPUTS
	- ordem: Set of feature which were selected.
	- maxcriterio: Value of criteria for the order calculated.
	"""
	pass


def selecaoVetorial(classes, K, metodo='exaustivo', criterio='J1'):
	""" Selects the best set of k features which separate the classes.

	INPUTS
	- classes: Python List of numpy matrixes of each class, row are patterns and columns are features
	- K: Number of features to be selected.
	- metodo: Type of method to be used: 'exaustivo', 'forward' or 'floating'
	- criterio: Criteria to be used to calculate the best set of features.

	OUTPUTS
	- ordem: Set of feature which were selected.
	- maxcriterio: Value of criteria for the order calculated.
	"""
	L = classes[0].shape[1]

	if K > L:
		raise AttributeError('Please, choose a smaller number of features to be selected.')

	if metodo.lower() == 'exaustivo':
		ordem, maxcriterio = exaustivosel(classes=classes, K=K, criterio=criterio)
	elif metodo.lower() == 'forward':
		ordem, maxcriterio = forwardsel(classes=classes, K=K, criterio=criterio)
	elif metodo.lower() == 'floating':
		ordem, maxcriterio = floatingsel(classes=classes, K=K, criterio=criterio)
	else:
		print('Método desconhecido!')
		return

	return ordem, maxcriterio
