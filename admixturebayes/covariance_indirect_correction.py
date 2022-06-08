import numpy as np
from covariance_simulation import Simulator
from copy import deepcopy
from scipy.stats import wishart, norm
from covariance_estimator import Estimator

class Sigma_proposal(object):
    
    def __init__(self, step_size=0.5,wait=6, threshold=0.4):
        self.step_size=step_size
        self.no_smallers=0
        self.reset_no_smallers=2
        self.wait=wait
        self.df=step_size
        self.threshold=threshold
        
    def is_wishart(self):
        if self.no_smallers>self.wait and (self.no_smallers%2)==((self.wait+1)%2):
            return True
        return False
    
    def is_wishart_era(self):
        return self.no_smallers>self.wait
    
    def __call__(self, emp_Sigma, implied_Sigma, Sigma):
        if self.is_wishart():
            print('Check wishart')
            df=emp_Sigma.shape[0]/self.df
            d_Sigma= wishart.rvs(df=df, scale=Sigma/df)-Sigma
            prop_Sigma=Sigma+np.absolute(d_Sigma)*np.sign(emp_Sigma-implied_Sigma)
            return prop_Sigma
            prop_Sigma=deepcopy(Sigma)
        else:
            prop_Sigma=Sigma+(emp_Sigma-implied_Sigma)*self.step_size
        return prop_Sigma
        
    def choose_next_and_adapt(self, old_distance, new_distance, old_Sigma, new_Sigma):
        if new_distance>old_distance:
            
            if self.is_wishart():
                #removedprin old_Sigma.shape
                self.df=min(0.5*self.df, 1)
            else:
                self.step_size*=0.8
                self.no_smallers+=1
            return old_distance, old_Sigma
        else:
            if self.is_wishart():
                self.no_smallers+=1
                self.df=min(5.0*self.df, 1)
                self.no_smallers=self.reset_no_smallers
            else:
                self.step_size*=1.1
            return new_distance, new_Sigma
        
        
class Chol_proposal(object):
    
    def __init__(self, step_size=0.5,wait=6, threshold=0.4):
        self.step_size=step_size
        self.no_smallers=0
        self.reset_no_smallers=2
        self.wait=wait
        self.df=step_size
        self.threshold=threshold
        
    def is_wishart(self):
        if self.no_smallers>self.wait and (self.no_smallers%2)==((self.wait+1)%2):
            return True
        return False
    
    def is_wishart_era(self):
        return self.no_smallers>self.wait
    
    
    def __call__(self, emp_Chol, implied_Chol, Chol):
        if self.is_wishart():
            print('Check wishart')
            print(Chol.shape)
            d_Chol= norm.rvs(size=Chol.shape)*self.df
            prop_Chol=Chol+np.absolute(d_Chol)*np.sign(emp_Chol-implied_Chol)
            return prop_Chol
            prop_chol=deepcopy(Sigma)
        else:
            prop_Chol=Chol+(emp_Chol-implied_Chol)*self.step_size
        return prop_Chol
        
    def choose_next_and_adapt(self, old_distance, new_distance, old_Sigma, new_Sigma, old_implied_Chol, new_implied_Chol):
        if new_distance>old_distance:
            
            if self.is_wishart():
                #removedprin old_Sigma.shape
                self.df=min(0.5*self.df, 1)
            else:
                self.step_size*=0.8
                self.no_smallers+=1
            return old_distance, old_Sigma, old_implied_Chol
        else:
            if self.is_wishart():
                self.no_smallers+=1
                self.df=min(5.0*self.df, 1)
                self.no_smallers=self.reset_no_smallers
            else:
                self.step_size*=1.1
            return new_distance, new_Sigma, new_implied_Chol
 
def turn_postive_definite(sym_mat):
    w,_ =np.linalg.eig(sym_mat)
    minw,maxw=np.min(w), np.max(w)
    if minw>0:
        add_on=0
    else:
        add_on=-minw+maxw/20.0
    print('add_on', add_on)
    print('w',w)
    print('sym_mat', sym_mat)
    print('new', sym_mat+add_on*np.identity(sym_mat.shape[0]))
    return sym_mat+add_on*np.identity(sym_mat.shape[0]), add_on
     
        
def status_print(i, proposal, emp_Sigma, implied_Sigma, Sigma, prop_Sigma, distances, prop_dist):
    to_print='Iteration '+str(i)+'\n'
    to_print+='Step size='+str(proposal.step_size)+'('+str(proposal.no_smallers)+")"+"\n"
    to_print+='Degrees of freedom='+str(emp_Sigma.shape[0]/proposal.df)+"\n"
    to_print+='{:<22} {:>8.5f}\n'.format('Old distance=', distances[-1])
    to_print+='{:<22} {:>8.5f}\n'.format('Proposed distance=', prop_dist)
    to_print+='emp Sigma and implied Sigma (diff= '+ str(np.linalg.norm(emp_Sigma-implied_Sigma))+"):"+'\n'
    to_print+=''.join(['{:<9.4f}'.format(s) for s in np.diag(emp_Sigma)])+' det='+str(np.linalg.det(emp_Sigma))+'\n'
    to_print+=''.join(['{:<9.4f}'.format(s) for s in np.diag(implied_Sigma)])+' det='+str(np.linalg.det(implied_Sigma))+'\n'
    to_print+='Diagonal elements of Sigma and proposed Sigma(diff= '+ str(np.linalg.norm(prop_Sigma-Sigma))+"):"+'\n'
    to_print+=''.join(['{:<9.4f}'.format(s) for s in np.diag(prop_Sigma)])+' det='+str(np.linalg.det(prop_Sigma))+'\n'
    to_print+=''.join(['{:<9.4f}'.format(s) for s in np.diag(Sigma)])+' det='+str(np.linalg.det(Sigma))+'\n'
    print(to_print)
    
def search_choleskys(xs,ns,no_its = 100, s=1, estimator=None, init_Sigma=None, Sim=None):
    #nss=np.tile(ns,s)
    if Sim is None:
        Sim=Simulator(nss, multiplier=s, estimator=estimator)
    no,N=xs.shape
    n=no-1
    emp_pijs=xs/ns
    emp_Sigma=estimator(xs,ns)#estimate_Sigma_wrapper(emp_pijs, reduce_method=reduce_method, method_of_weighing_alleles=method_of_weighing_alleles)
    emp_Sigma, a=turn_postive_definite(emp_Sigma)
    emp_Chol=np.linalg.cholesky(emp_Sigma)
    if init_Sigma is None:
        Chol=emp_Chol
    else:
        Chol=np.linalg.cholesky(init_Sigma)
    implied_Sigma=Sim.get_implied_Sigma_from_chol(Chol)+np.identity(emp_Chol.shape[0])*a
    implied_Chol=np.linalg.cholesky(implied_Sigma)
    distances=[np.linalg.norm(emp_Sigma-implied_Sigma)]
    Proposal=Chol_proposal()
    for i in range(no_its):
        prop_Chol=Proposal(emp_Chol, implied_Chol, Chol)
        try:
            implied_Sigma=Sim.get_implied_Sigma_from_chol(prop_Chol)+np.identity(prop_Chol.shape[0])*a
            prop_dist=np.linalg.norm(emp_Sigma-implied_Sigma)
            new_implied_Chol=np.linalg.cholesky(implied_Sigma)
        except np.linalg.linalg.LinAlgError:
            prop_dist=float('inf')
            new_implied_Chol=implied_Chol
        
        status_print(i, Proposal, emp_Sigma, implied_Sigma, Chol.dot(Chol.T), prop_Chol.dot(prop_Chol.T), distances, prop_dist)
        new_dist, Chol, implied_Chol= Proposal.choose_next_and_adapt(distances[-1], prop_dist, Chol,prop_Chol, implied_Chol, new_implied_Chol)
        distances.append(new_dist)

    return Chol.dot(Chol.T), distances[-1]


def search_sigmas(xs,ns,no_its = 100, s=1, estimator=None, init_Sigma=None, Sim=None):
    #nss=np.tile(ns,s)
    if Sim is None:
        Sim=Simulator(nss, multiplier=s, estimator=estimator)
    no,N=xs.shape
    n=no-1
    emp_pijs=xs/ns
    emp_Sigma=estimator(xs,ns)#estimate_Sigma_wrapper(emp_pijs, reduce_method=reduce_method, method_of_weighing_alleles=method_of_weighing_alleles)
    if init_Sigma is None:
        Sigma=emp_Sigma
    else:
        Sigma=init_Sigma
    implied_Sigma=Sim.get_implied_Sigma(Sigma)
    distances=[np.linalg.norm(emp_Sigma-implied_Sigma)]
    Proposal=Sigma_proposal()
    for i in range(no_its):
        prop_Sigma=Proposal(emp_Sigma, implied_Sigma, Sigma)
        try:
            implied_Sigma=Sim.get_implied_Sigma(prop_Sigma)
            prop_dist=np.linalg.norm(emp_Sigma-implied_Sigma)
        except np.linalg.linalg.LinAlgError:
            prop_dist=float('inf')
        
        status_print(i, Proposal, emp_Sigma, implied_Sigma, Sigma, prop_Sigma, distances, prop_dist)
        new_dist, new_Sigma= Proposal.choose_next_and_adapt(distances[-1], prop_dist, Sigma,prop_Sigma)
        distances.append(new_dist)
        Sigma=new_Sigma

    return Sigma, distances[-1]
    
class IndirectEstimator(Estimator):
    

    def __init__(self, 
                 no_its=100,
                 s=1,
                 initial_Sigma_generator=None,
                 estimator=None,
                 simulator=None,
                 ):
        super(IndirectEstimator, self).__init__(reduce_also=estimator.reduce_also)
        self.s=s
        self.no_its=no_its
        self.initial_Sigma_generator=initial_Sigma_generator
        self.initialize_Sigma()
        
        self.estimator=estimator
        self.simulator=simulator
    
    
        
    def __call__(self, xs,ns, extra_info={}):
        cov,val= search_choleskys(xs,
                             ns,
                             no_its = self.no_its, 
                             s=self.s, 
                             estimator=self.estimator, 
                             init_Sigma=self.initial_Sigma, 
                             Sim=self.simulator)
        self.fitted_value=val
        return cov
