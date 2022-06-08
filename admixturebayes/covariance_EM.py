import numpy as np
from covariance_estimator import Estimator


def update_pop(Sigma, hidden_states, i, freqs, vars, p0s):
    Sigma_11=Sigma[i,i]
    Sigma_21=np.delete(Sigma[:,i],i,0)
    #removedprin Sigma_21
    Sigma_12=Sigma_21.T
    Sigma_22i=np.linalg.inv(np.delete(np.delete(Sigma, i,0),i,1))
    hidden_states_small=np.delete(hidden_states,i,0)
    #removedprin hidden_states_small[:,1]
    N=hidden_states.shape[1]
    Sigma_prod=np.dot(Sigma_12, Sigma_22i)
    j=1
    #removedprin p0s[j]+np.dot(Sigma_prod , hidden_states_small[:,j] - p0s[j]  )
    means=np.array([p0s[j]+np.dot(Sigma_prod , hidden_states_small[:,j] - p0s[j]  ) for j in range(N)])
    #removedprin means
    var_base=Sigma_11+np.dot(Sigma_prod, Sigma_21.T)
    #removedprin var_base
    vars_normal=np.array([var_base*p0*(1-p0) for p0 in p0s])
    means2=freqs[i,:].flatten()
    vars2=vars[i,:].flatten()
    #removedprin means.shape, vars_normal.shape, means2.shape, vars2.shape
    return (means/vars_normal+means2/vars2)/(1.0/vars_normal+1.0/vars2)

def update_sigma(hidden_states, p0s):
    scaled=(hidden_states-p0s.T)/np.sqrt(p0s.T*(1.0-p0s.T))
    return np.dot(scaled,scaled.T)/hidden_states.shape[1]

def heuristic_p0(xs,ns, alpha=1.0):
    freqs=xs/ns
    freqs=np.clip(freqs,0.01,0.99)
    return freqs[0,:]*alpha+np.mean(freqs[1:], axis=0)*(1-alpha)

def initor(a):
    if not isinstance(a, str):
        return a[0]
    else:
        return a

def em(xs, ns, initial_Sigma=None, p0s=None, alpha=1.0, maxiter=100):
    freqs=xs[1:,:]/ns[1:,:]
    freqs=np.clip(freqs,0.01,0.99)
    print(np.min(freqs))
    hidden_states=xs[1:,:]/ns[1:,:]
    hidden_states=np.clip(hidden_states,0.01,0.99)
    if p0s is None:
        p0s=heuristic_p0(xs, ns, alpha=alpha)
    p0s=p0s.clip(0.01,0.99)
    if initial_Sigma is None:
        Sigma=update_sigma(hidden_states, p0s)
    else:
        Sigma=initial_Sigma
    vars=freqs*(1.0-freqs)
    n=xs.shape[0]-1
    count=0
    old_Sigma=Sigma
    first_dist=0
    while count<maxiter:
        for i in range(n):
            Sigma=update_sigma(hidden_states, p0s)
            hidden_states[i,:]=update_pop(Sigma, hidden_states, i, freqs, vars, p0s).T
        dist=np.linalg.norm(old_Sigma-Sigma)
        print('{}: {:.6f}'.format(count,dist))
        old_Sigma=Sigma
        if dist<1e-5:
            break
            
        
        count+=1
        
    return Sigma

class EmEstimator(Estimator):
    

    def __init__(self, 
                 maxiter=100,
                 alpha=1.0,
                 initial_Sigma_generator=None,
                 ):
        super(EmEstimator, self).__init__(reduce_also=True)
        self.alpha=alpha
        self.maxiter=maxiter
        self.initial_Sigma_generator=initial_Sigma_generator
        self.initialize_Sigma()

        
    def __call__(self, xs,ns):
        return em(xs, ns, initial_Sigma=self.initial_Sigma, alpha=self.alpha, maxiter=self.maxiter)
