from estimators import Estimator, initor
from numpy import array, mean, zeros, diag, sum, savetxt,nan, isnan, nanmean
from numpy import sum as npsum
import warnings

from numpy import insert, identity

def reduce_covariance(covmat, subtracted_population_index):
    reducer=insert(identity(covmat.shape[0]-1), subtracted_population_index, -1, axis=1)
    return reducer.dot(covmat).dot(reducer.T)

default_scale_dic={'None':'None',
                   'Jade-o':'outgroup_sum', 
                   'Jade':'average_sum',
                   'outgroup_sum':'outgroup_sum',
                   'outgroup_product':'outgroup_product',
                   'average_sum':'average_sum',
                   'average_product':'average_product'}

def nan_divide(dividend, divisor):
    res=zeros(dividend.shape)
    for i in range(dividend.shape[0]):
        for j in range(dividend.shape[1]):
            if divisor[i,j]==0:
                res[i,j]=nan
            else:
                res[i,j]=dividend[i,j]/divisor[i,j]
    return res

def nan_inner_product(a,b):
    N=len(a)
    nans=0
    res=0
    for ai,bi in zip(a,b): 
        if isnan(ai) or isnan(bi):
            nans+=1
        else:
            res+=ai*bi
    if N==nans:
        warnings.warn('There is an entry in the covariance matrix that is set to 0 because all the relevant data was nan.', UserWarning)
        return 0
        
    return res/(N-nans)

def nan_product(A,B):
    res=zeros((A.shape[0], B.shape[1]))
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            res[i,j]=nan_inner_product(A[i,:], B[:,j])
    return res

def m_scaler(allele_freqs):
    s=nanmean(allele_freqs, axis=0)
    scaler=nanmean(s*(1.0-s))
    return scaler

def avg_var(ps):
    return sum((p*(1-p) for p in ps))/float(len(ps))

def heterogeneity(allele_frequency, pop_size, type_of_scaling='unbiased'):
    if type_of_scaling=='mle':
        mult=2.0/float(len(allele_frequency))#*float(pop_size)/float(pop_size-1)
    else:
        mult=2.0/float(len(allele_frequency))*float(pop_size)/float(pop_size-1)
    return sum([p*(1-p) for p in allele_frequency])*mult

def B(allele_frequency, pop_size, type_of_scaling='unbiased'):
    return heterogeneity(allele_frequency, pop_size, type_of_scaling)/2.0/float(pop_size)

def var(p,n, type_of_scaling='unbiased'):
    if type_of_scaling=='mle':
        subtract=0
    else:
        subtract=1
    entries=array([pi*(1-pi)/(ni-subtract) for pi,ni in zip(p,n) if ni>subtract])
    return mean(entries)

def reduced_covariance_bias_correction(p,n,n_outgroup=0, type_of_scaling='unbiased'):
    Bs=[]
    for pi,ni in zip(p,n):
        Bs.append(var(pi,ni, type_of_scaling))
    outgroup_b=Bs.pop(n_outgroup)
    return diag(array(Bs))+outgroup_b

class ScaledEstimator(Estimator):
    
    def __init__(self,
                 reduce_method=['outgroup','average','None'],
                 scaling=['None','outgroup_sum', 'outgroup_product', 'average_outgroup', 'average_product','Jade','Jade-o'],
                 variance_correction=['None','unbiased','mle'],
                 add_variance_correction_to_graph=False,
                 prefix_for_saving_variance_correction='',
                 save_variance_correction=True,
                 nodes=None):
        super(ScaledEstimator, self).__init__(reduce_also=True)
        self.scaling=initor(scaling)
        self.variance_correction=initor(variance_correction)
        self.add_variance_correction_to_graph=add_variance_correction_to_graph
        self.prefix_for_saving_variance_correction=prefix_for_saving_variance_correction
        self.nodes=nodes
        self.save_variance_correction=save_variance_correction
        self.bias_c_weight=default_scale_dic[scaling]
        self.reduce_method=reduce_method
        
    def subtract_ancestral_and_get_outgroup(self,p):
        n_outgroup=self.get_reduce_index()
        return p-p[n_outgroup,:], n_outgroup
        
    def __call__(self, xs, ns, extra_info={}):
        if 0 in ns:
            warnings.warn('There were 0s in the allele-totals, inducing nans and slower estimation.', UserWarning)
            ps=nan_divide(xs, ns)
        else:
            ps=xs/ns 
        return self.estimate_from_p(ps, ns=ns, extra_info=extra_info)
    
        
    def estimate_from_p(self, p, ns=None, extra_info={}):
        
        p2,n_outgroup = self.subtract_ancestral_and_get_outgroup(p)
        
        if npsum(isnan(p2))>0:
            warnings.warn('Nans found in the allele frequency differences matrix => slower execution', UserWarning)
            m=nan_product(p2, p2.T)
        else:
            m=p2.dot(p2.T)/p2.shape[1]
        assert m.shape[0]<1000, 'sanity check failed, because of wrongly transposed matrices'
        
        scaling_factor=m_scaler( p)
        extra_info['m_scale']=scaling_factor
        m=m/scaling_factor
        if self.reduce_also:
            m=reduce_covariance(m, n_outgroup)
            if self.variance_correction!='None':
                assert ns is not None, 'Variance correction needs a ns-matrix specified'
                b=reduced_covariance_bias_correction(p, ns, n_outgroup, type_of_scaling=self.variance_correction)/scaling_factor
                if self.add_variance_correction_to_graph:
                    if self.save_variance_correction:
                        savetxt(self.prefix_for_saving_variance_correction+'variance_correction.txt', b)
                else:
                    m=m-b
        return m