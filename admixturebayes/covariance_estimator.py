#from numpy import array, mean, zeros, diag, sum, arcsin, sqrt

def initor(a):
    if not isinstance(a, str):
        return a[0]
    else:
        return a
    
class Estimator(object):
    
    def __init__(self,
                 reduce_also=True,
                 arcsin_transform=False):
        self.arcsin_transform=arcsin_transform
        self.reduce_also=reduce_also
        self.fitted_value=None
        self.initial_Sigma_generator=None
        self.initial_Sigma=None
        
    def initialize_Sigma(self):
        if self.initial_Sigma_generator is None:
            self.initial_Sigma=None
        else:
            self.initial_Sigma=self.initial_Sigma_generator()  
        
    def get_fitted_value(self):
        '''
        Some initia
        '''
        return self.fitted_value
        
    def get_reduce_index(self):
        return 0
#         n_outgroup=next((n for n, e in enumerate(self.full_nodes) if e==self.outgroup))
#         return n_outgroup
    
    def __call__(self,xs,ns, names=None, extra_info={}):
        '''
        Should return the covariance estimate and possibly set the variable fitted_value
        '''
        pass
    
class RepeatEstimator(Estimator):
    
    def __init__(self, est, reps=100):
        
        super(RepeatEstimator, self).__init__(reduce_also=est.reduce_also)
        self.est=est
        self.reps=reps
        
    def __call__(self, xs,ns, extra_info={}):
        vals=[]
        covs=[]
        for i in range(self.reps):
            self.est.initialize_Sigma()
            cov=self.est(xs,ns, extra_info=extra_info)
            vals.append(self.est.get_fitted_value())
            covs.append(cov)
        
        min_val=vals[0]
        min_index=0
        for i,(cov, val) in enumerate(zip(covs, vals)):
            print(val, ':', cov)
            if val is not None and val < min_val:
                min_index=i
                min_val=val
        self.fitted_value=vals[min_index]
        return covs[min_index]
            
        
