import numpy as np
from scipy.signal import welch
from scipy.stats import skew, kurtosis
from scipy.interpolate import Rbf
from itertools import permutations, combinations
import matplotlib.pyplot as plt 
import spectrum

def featureExtraction(signal, Fs, field='sleep', bands=0, order=2):
	"""Feature extraction function based on pertinent features for sleep or hrv analysis
	
	INPUTS
	- signal: Matrix N by M, with N signals and M samples.
	- Fs: Sample rate of the signals.
	- field: What kind of features will be calculated: "sleep" or "hrv".
	- bands: Number of eeg bands to be calculated. Only for sleep analysis.

	OUTPUTS
	- data: matrix of features in each column.
	
	"""

	N, M = signal.shape # N patterns and M samples
	if field.lower() == 'sleep': # Calculates relevant features for sleep analysis
		data = featureExtractionSleep(signal, Fs, bands=0)
	elif field.lower() == 'hrv': # Calculates relevant features for HRV analysis
		data = featureExtractionHrv(signal, Fs, order=2)

	return data

def featureExtractionSleep(signal, Fs, bands=0):
	"""Feature extraction function based on pertinent features for sleep or hrv analysis
	
	INPUTS
	- signal: Matrix N by M, with N signals and M samples.
	- Fs: Sample rate of the signals.
	- field: What kind of features will be calculated: "sleep" or "hrv".
	- bands: Number of eeg bands to be calculated. Only for sleep analysis.

	OUTPUTS
	- data: matrix of features in each column.
	
	"""

	N, M = signal.shape # N patterns and M samples

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
	data[7, :] = data[4, :]/data[2, :]                                 # Mobilidade Estatistica 
	data[8, :] = np.sqrt(data[6, :]/data[4, :] - data[7, :])           # Complexidade Estatistica 
	data[9, :] = skew(signal, axis=1)                                  # Obliquidade
	data[10, :] = kurtosis(signal, axis=1)                             # Kurtose
	for n in range(0, N):	
		# Caracteristicas Espectrais
		f, pxx = welch(signal[n, :], fs=Fs, window='hanning', noverlap=50)

		data[11, n] = np.sum(f*pxx)/np.sum(pxx)                                    # Frequencia Central
		data[12, n] = pxx[np.where(f-data[11, :] == min(f-data[11, :]))[0]]        # Potencia na frequencia central
		data[13, n] = np.sum((f-data[11, n])*pxx)/np.sum(pxx)                      # Largura de Banda
		data[14, n] = f[int(min(f[np.cumsum(pxx) > np.sum(pxx)*0.9]))]             # Frequencia de margem
	"""
	for c in range(13, L):
		data[c, n] = 
	"""	
	
	return data.T

def featureExtractionHrv(signal, Fs, order=2):
	"""Feature extraction function based on pertinent features for sleep or hrv analysis
	
	INPUTS
	- signal: Matrix N by M, with N signals and M samples.
	- Fs: Sample rate of the signals.
	- field: What kind of features will be calculated: "sleep" or "hrv".
	- bands: Number of eeg bands to be calculated. Only for sleep analysis.

	OUTPUTS
	- data: matrix of features in each column.
	
	"""

	N, M = signal.shape # N patterns and M samples
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

		# Caracteristicas Espectrais
		
		pxx = pburg(signal[n, :], order=order, sampling=Fi)
		freq = np.linspace(0, Fi, pxx.shape)

		data[4, n] = np.sum(freq*pxx)/np.sum(pxx)                                    # Frequencia Central
		data[5, n] = pxx[np.where(freq-data[11, :] == min(freq-data[11, :]))[0]]        # Potencia na frequencia central
		data[6, n] = np.sum((freq-data[11, n])*pxx)/np.sum(pxx)                      # Largura de Banda
		data[7, n] = freq[int(min(freq[np.cumsum(pxx) > np.sum(pxx)*0.9]))]             # Frequencia de margem

	return data.T