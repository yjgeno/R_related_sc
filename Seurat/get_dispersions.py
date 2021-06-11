import numpy as np
import matplotlib.pyplot as plt


def get_Mean_Var(X):
    if isinstance(X, np.ndarray):
        mean = np.mean(X, axis=0, dtype=np.float64)
        mean_sq = np.multiply(X, X).mean(axis=0, dtype=np.float64)
        var = mean_sq - mean ** 2
        var *= X.shape[0] / (X.shape[0] - 1) #unbiased variance
    return mean, var
  
  
def get_Dispersions(X, log = True, plot = False): #X = ada.X, log space X
    X = np.expm1(X) #to original counts
    mean, var = get_Mean_Var(X)

    mean[mean == 0] = 1e-12  
    dispersion = var / mean
    dispersion[dispersion == 0] = np.nan #for log   
    CV = dispersion / np.sqrt(var)
    
    if log:
        mean = np.log2(mean)
        CV= np.log2(CV)
        dispersion= np.log2(dispersion)
         
    if plot:
        plt.figure(figsize=(12,6))
        plt.subplot(1, 2, 1)
        plt.scatter(mean, CV, marker='o', s=3)
        plt.xlabel('(log) mean')
        plt.ylabel('(log) CV')
        
        plt.subplot(1, 2, 2)
        plt.scatter(mean, dispersion, marker='o', s=3)
        plt.xlabel('(log) mean')
        plt.ylabel('(log) dispersion')
        plt.show()
    
    return dispersion, CV 
  
  
#test
#dispersion, CV = get_Dispersions(ada.X, log = True, plot=True) 

  
  
  
  
  
