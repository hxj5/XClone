# xmath.py

import logging
import numpy as np
import scipy as sp
from scipy.special import logsumexp


def cal_log_lik(emm_prob_log, posterior_mtx_log):
    """
    Function:

    INPUT:
    -------
    c * g * states
    c * g * states
    c can be celltype based | can be cell based
    """
    emm_dim = emm_prob_log.ndim
    pos_dim = posterior_mtx_log.ndim

    if emm_dim != 3 | pos_dim != 3:
        logging.warning("dim mismatch: emm_dim = %d, pos_dim = %d." % 
            (emm_dim, pos_dim))
        return None

    log_lik = logsumexp(emm_prob_log + posterior_mtx_log, axis = 2).sum(axis = 1).sum()
    return log_lik


def loglik_amplify(X, axis = -1):
    """
    Amplify the log likelihood matrix by subtract the maximum.
    X should be numpy.array, otherwise will be transformed to it.
    
    Example
    -------
    X = np.random.rand(3, 5, 8)
    loglik_amplify(X, axis = 1)
    """    
    if type(X) == np.matrix:
        X = X.A
    
    X_max = np.max(X, axis = axis, keepdims = True)
    return X - X_max


def normalize(X, axis = -1):
    """
    Normalization of tensor with sum to 1.
    X should be numpy.array, otherwise will be transformed to it.
    
    Example
    -------
    X = np.random.rand(3, 5, 8)
    normalize(X, axis = 1)
    """    
    if type(X) == np.matrix:
        X = X.A
    if sp.sparse.issparse(X):
        X = X.A
    
    X_sum = np.sum(X, axis = axis, keepdims = True)
    return X / (X_sum + 1e-8)


# refer to https://stackoverflow.com/questions/21610198/runtimewarning-divide-by-zero-encountered-in-log
def safe_log(x, eps = 1e-6, *args, **kwargs):
    assert 0 < eps < 1
    #result = np.where(x > eps, x, np.log(eps))     
    #np.log(result, out = result, where = result > 0)
    x_tmp = np.where(x > eps, x, eps)
    result = np.log(x_tmp)
    return result

