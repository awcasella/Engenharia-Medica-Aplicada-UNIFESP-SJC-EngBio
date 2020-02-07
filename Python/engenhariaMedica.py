import numpy as np
from scipy.signal import welch
from scipy.stats import skew, kurtosis
from scipy.interpolate import Rbf
from itertools import permutations, combinations
import matplotlib.pyplot as plt 

def featureExtraction(signal, Fs, method='welch', field='sleep', bands=0):
	"""Feature extraction function based on pertinent features for sleep or hrv analysis
	
	INPUTS
	- signal: Matrix N by M, with N signals and M samples.
	- Fs: Sample rate of the signals.
	- method: Spectral method to be used. "welch" or "burg"
	- field: What kind of features will be calculated: "sleep" or "hrv".
	- bands: Number of eeg bands to be calculated. Only for sleep analysis.

	OUTPUTS
	- data: matrix of features in each column.
	
	"""

	N, M = signal.shape
	if field.lower() == 'sleep': # Calculates relevant features for sleep analysis
		#L = 15 + bands.shape[0]
		L = 15
		data = np.zeros([L, N])
		
		# Caracteristicas temporais
		data[0, :] = np.average(signal, axis=1)                            # Media
		data[1, :] = np.average(abs(signal), axis=1)                       # Media retificada
		data[2, :] = np.var(signal, axis=1, ddof=1)                        # Variancia
		data[3, :] = np.average(np.diff(signal, n=1, axis=1), axis=1)      # Media da primeira derivada
		data[4, :] = np.var(np.diff(signal, n=1, axis=1), axis=1, ddof=2)  # Variancia da primeira derivada
		data[5, :] = np.average(np.diff(signal, n=2, axis=1), axis=1)      # Media da segunda derivada
		data[6, :] = np.var(np.diff(signal, n=2, axis=1), axis=1, ddof=2)  # Variancia da segunda derivada
		data[7, :] = data[4, :]/data[2, :]                               # Mobilidade Estatistica 
		data[8, :] = np.sqrt(data[6, :]/data[4, :] - data[7, :])        # Complexidade Estatistica 
		data[9, :] = skew(signal, axis=1)                                  # Obliquidade
		data[10, :] = kurtosis(signal, axis=1)                             # Kurtose
		for n in range(0, N):	
			# Caracteristicas Espectrais
			if method.lower() == 'welch':
				f, pxx = welch(signal[n, :], fs=Fs, window='hanning', noverlap=50)
			elif method.lower() == 'burg':
				print('Em construcao... Usando metodo welch')
				f, pxx = welch(signal[n, :], fs=Fs, window='hanning', noverlap=50, axis=1)

			data[11, n] = np.sum(f*pxx)/np.sum(pxx)                                    # Frequencia Central
			data[12, n] = pxx[np.where(f-data[11, :] == min(f-data[11, :]))[0]]        # Potencia na frequencia central
			data[13, n] = np.sum((f-data[11, n])*pxx)/np.sum(pxx)                      # Largura de Banda
			data[14, n] = f[int(min(f[np.cumsum(pxx) > np.sum(pxx)*0.9]))]             # Frequencia de margem
		"""
		for c in range(13, L):
			data[c, n] = 
		"""	
	elif field.lower() == 'hrv': # Calculates relevant features for HRV analysis
		Fi = Fs # interpolation sample rate
		x = np.linspace(1, M, M)
		xx = np.linspace(1, M, M*Fi)
		L = 4

		data = np.zeros([L, N])
		for n in range(0, N):
			f = Rbf(x, sinal[n, :], function='cubic')
			sinalInterpolado = f(xx)
			data[0, n] = sum((60/(sinalInterpolado)).T).T                    # Heart rate
			data[1, n] = np.std(sinalInterpolado)                            # SDNN
			data[2, n] = np.average(np.diff(sinalInterpolado, n=1)**2)**0.5  # rMSSD
			data[3, n] = sum((np.diff(sinalInterpolado, n=1)).T).T           # PNN50

	return data.T

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


"""
AULA 28
"""

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
	# FDA 
	# 
	# INPUTS:
	# - classes: Celula de classes.
	# - n: Numero de caracteristicas do novo espaco.
	# 
	# OUTPUTS:
	# - Y: Novo espaço otimizado
	L = classes[0].shape[0]
	Sw, Sb = scatter(classes, L)                   # Ayrton, acho que aqui eh L.
	eigval, eigvec = np.linalg.eig(np.linalg.inv(Sw).dot(Sb))

	eigvec, eigval = autoordenar(eigvec, eigval)

	C = len(classes)                               # Numero de classes.
	A = eigvec[:, 0:C-1].T

	X = classes[0]

	for n in range(1, C):
		X = np.concatenate((X, classes[n]), axis=1)


	Y = A.dot(X)

	return Y
