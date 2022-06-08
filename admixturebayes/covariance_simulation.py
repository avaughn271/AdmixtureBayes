from scipy.stats import norm, uniform, binom, multivariate_normal
import numpy as np

class Simulator(object):
    
    def __init__(self, ns,multiplier=1, estimator=None, fixed_seed=True,  load_from_file='', locus_filter=None):
        self.ns=np.tile(ns, multiplier)
        self.n=self.ns.shape[0]-1
        self.N=self.ns.shape[1]
        print('self.N=', self.N)
        print('self.n=', self.n)
        self.fixed_seed=fixed_seed
        self.initialize_sims(load_from_file)
        self.initialize_nvals()
        self.estimator=estimator
        self.locus_filter=locus_filter
        
    def initialize_nvals(self):
        self.nvals=np.unique(self.ns)
        self.ids=[]
        for val in self.nvals:
            self.ids.append(np.where(self.ns==val))
            
        
        
    def get_xs(self, Sigma):
        
        L=np.linalg.cholesky(Sigma)
        return self.get_xs_from_chol(L)
    
    def get_xs_from_chol(self, L):
        pijs=np.dot(L, self.Us)+self.p0s
        trunc_pijs=np.clip(pijs,0,1)
        trunc_pijs=np.insert(trunc_pijs,0,self.p0s,axis=0)
        x_ijs=self.qbinom(trunc_pijs)
        x_ijs=adjust_scipy_error(trunc_pijs, self.ns, x_ijs)

        return x_ijs
    
    def get_implied_Sigma(self, Sigma):
        if not self.fixed_seed:
            self.initialize_sims(False)
        x_ij=self.get_xs(Sigma)
        print(x_ij.shape)
        if self.locus_filter is not None:
            xs,ns=self.locus_filter.apply_filter(x_ij/self.ns,self.ns)
        print(xs.shape)
        cov=self.estimator(xs*ns, ns) #if this line fails, it could be because estimator has default value None
        return cov
    
    def get_implied_Sigma_from_chol(self, Chol):
        if not self.fixed_seed:
            self.initialize_sims(False)
        x_ij=self.get_xs_from_chol(Chol)
        if self.locus_filter is not None:
            xs,ns=self.locus_filter.apply_filter(x_ij/self.ns,self.ns)
        print(xs.shape)
        cov=self.estimator(xs*ns, ns)
        return cov
        
            
    def initialize_sims(self, load_from_file):
        if load_from_file:
            self.initialize_from_file(filename_prefix=load_from_file)
            return
        print('SIMULATING')
        self.p0s=uniform.rvs(size=self.N)
        #self.x0s=self.qbinom(self.p0s)
        self.Us=norm.rvs(size=(self.n,self.N))
        self.Us*=np.sqrt(self.p0s*(1.0-self.p0s))
        self.Vs=uniform.rvs(size=(self.n+1,self.N))
        
    def initialize_from_file(self, filename_prefix):
        self.ns= np.loadtxt(filename_prefix+'ns.txt')
        self.p0s= np.loadtxt(filename_prefix+'p0s.txt')
        self.Us= np.loadtxt(filename_prefix+'Us.txt')
        self.Vs= np.loadtxt(filename_prefix+'Vs.txt')
        self.n=self.ns.shape[0]-1
        self.N=self.ns.shape[1]
        
    def qbinom(self, trunc_pijs):
        res=np.zeros(trunc_pijs.shape)
        for i,n in enumerate(self.nvals):
                res[self.ids[i]]=binom.ppf(self.Vs[self.ids[i]], n=n, p=trunc_pijs[self.ids[i]])
        return res
    
    def save_to_file(self, filename_prefix):
        np.savetxt(filename_prefix+'ns.txt', self.ns)
        np.savetxt(filename_prefix+'p0s.txt', self.p0s)
        np.savetxt(filename_prefix+'Us.txt',self.Us)
        np.savetxt(filename_prefix+'Vs.txt',self.Vs)        
    
def adjust_scipy_error(pijs, ns, xs):
    xs[np.where(pijs<1e-10)]=0
    ids=np.where(pijs>1-1e-10)
    xs[ids]=ns[ids]
    return xs
    