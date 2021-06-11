import numpy as np
import matplotlib.pyplot as plt


def get_Mean_Var(X):
    if isinstance(X, np.ndarray):
        mean = np.mean(X, axis=0, dtype=np.float64)
        mean_sq = np.multiply(X, X).mean(axis=0, dtype=np.float64)
        var = mean_sq - mean ** 2
        var *= X.shape[0] / (X.shape[0] - 1) #unbiased variance
    return mean, var
  
  
def get_Dispersions(X, plot = False):
    #X = ada.X #log space X
    X = np.expm1(X)
    mean, var = get_Mean_Var(X)

    mean[mean == 0] = 1e-12  
    dispersion = var / mean
    dispersion[dispersion == 0] = np.nan #for log   
    CV = dispersion / np.sqrt(var)
    
    if plot:
        plt.scatter(np.log2(mean), np.log2(CV), marker='o', s=3)
        plt.xlabel('log2 mean')
        plt.ylabel('log2 CV')
        #plt.scatter(np.log2(mean), np.log2(var), marker='o')
        plt.show()
    
    return dispersion, CV #np.log(dispersion)
  
  
#test
#get_Dispersions(ada.X)[:5]
  
  
  
  
  
