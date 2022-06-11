from covariance_estimator import Estimator, initor
from numpy import array, mean, zeros, diag, sum, sqrt, savetxt,nan, isnan, nanmean
from numpy import sum as npsum
from reduce_covariance import reduce_covariance
import warnings

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

def m_scaler(scale_type, allele_freqs, n_outgroup=None):
    if scale_type=='None' or scale_type=='Jade' or scale_type=='Jade-o':
        return 1.0
    if scale_type.startswith('outgroup'):
        s=allele_freqs[n_outgroup,:]
    elif scale_type.startswith('average'):
        s=nanmean(allele_freqs, axis=0)
    else:
        scaler=1.0
    if scale_type.endswith('product'):
        mu=nanmean(s)
        scaler=mu*(1.0-mu)
    elif scale_type.endswith('sum'):
        scaler=nanmean(s*(1.0-s))
    #removedprin 'm_scale', scaler
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

def adjuster(Bs):
    m=len(Bs)
    res=diag(array(Bs))
    res=res-array(Bs)/m
    res=(res.T-array(Bs)/m).T
    res=res+sum(Bs)/m**2
    return res

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
    
    
    

def bias_correction(m, p, pop_sizes, n_outgroup=None, type_of_scaling='unbiased'):
    Bs=[B(prow, pop_size, type_of_scaling=type_of_scaling) for prow, pop_size in zip(p, pop_sizes)]
    adjusting_matrix=adjuster(Bs)
    return adjusting_matrix

class ScaledEstimator(Estimator):
    
    def __init__(self,
                 reduce_method=['outgroup','average','None'],
                 scaling=['None','outgroup_sum', 'outgroup_product', 'average_outgroup', 'average_product','Jade','Jade-o'],
                 reduce_also=True,
                 variance_correction=['None','unbiased','mle'],
                 jade_cutoff=1e-5,
                 bias_c_weight='default',
                 add_variance_correction_to_graph=False,
                 prefix_for_saving_variance_correction='',
                 save_variance_correction=True,
                 nodes=None):
        super(ScaledEstimator, self).__init__(reduce_also=reduce_also)
        self.scaling=initor(scaling)
        self.variance_correction=initor(variance_correction)
        self.jade_cutoff=jade_cutoff
        self.add_variance_correction_to_graph=add_variance_correction_to_graph
        self.prefix_for_saving_variance_correction=prefix_for_saving_variance_correction
        self.nodes=nodes
        self.save_variance_correction=save_variance_correction
        if bias_c_weight=='default':
            self.bias_c_weight=default_scale_dic[scaling]
        else:
            self.bias_c_weight=bias_c_weight
        self.reduce_method=reduce_method
        
    def subtract_ancestral_and_get_outgroup(self,p):
        if self.reduce_method=='outgroup':
            n_outgroup=self.get_reduce_index()
            #removedprin n_outgroup
            return p-p[n_outgroup,:], n_outgroup
        elif self.reduce_method=='average':
            n_outgroup=self.get_reduce_index()
            total_mean2=nanmean(p, axis=0)
            return p-total_mean2, n_outgroup
        else:
            return p, None
        
    def __call__(self, xs, ns, extra_info={}):
        if 0 in ns:
            warnings.warn('There were 0s in the allele-totals, inducing nans and slower estimation.', UserWarning)
            ps=nan_divide(xs, ns)
        else:
            ps=xs/ns 
        return self.estimate_from_p(ps, ns=ns, extra_info=extra_info)
    
        
    def estimate_from_p(self, p, ns=None, extra_info={}):
        #p=reorder_rows(p, self.names, self.full_nodes)
        
        p2,n_outgroup = self.subtract_ancestral_and_get_outgroup(p)
        
        
    
        if self.scaling=='Jade':
            mu=mean(p, axis=0)
            
            i=array([v > self.jade_cutoff and v<1.0-self.jade_cutoff for v in mu ])
            p2=p2[:,i]/sqrt(mu[i]*(1.0-mu[i]))
        elif self.scaling=='Jade-o':
            mu=p[n_outgroup,:]
            
            i=array([v > self.jade_cutoff and v<1.0-self.jade_cutoff for v in mu ])
            #p=p[:,i]
            p2=p2[:,i]/sqrt(mu[i]*(1.0-mu[i]))
        if npsum(isnan(p2))>0:
            warnings.warn('Nans found in the allele frequency differences matrix => slower execution', UserWarning)
            m=nan_product(p2, p2.T)
        else:
            m=p2.dot(p2.T)/p2.shape[1]
        assert m.shape[0]<1000, 'sanity check failed, because of wrongly transposed matrices'
        
        
        scaling_factor=m_scaler(self.scaling, p, n_outgroup)
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
        elif self.variance_correction!='None':
            m=m/m_scaler(self.scaling, p, n_outgroup)    
            extra_info['m_scale']=m_scaler(self.scaling, p, n_outgroup)     
            if ns is None:
                warnings.warn('No variance reduction performed due to no specified sample sizes', UserWarning)
            elif isinstance(ns, int):
                pop_sizes=[ns]*p2.shape[0]
                changer=bias_correction(m,p, pop_sizes,n_outgroup, type_of_scaling=self.variance_correction)/m_scaler(self.bias_c_weight, p, n_outgroup)
                #removedprin 'm',reduce_covariance(m,n_outgroup)
                #removedprin 'changer', reduce_covariance(changer, n_outgroup)
                if self.add_variance_correction_to_graph:
                    m-=changer
                if self.save_variance_correction:
                    savetxt(self.prefix_for_saving_variance_correction+'variance_correction.txt', changer)
            else:
                warnings.warn('assuming the same population size for all SNPs', UserWarning)
                pop_sizes=mean(ns, axis=1)
                changer=bias_correction(m,p, pop_sizes,n_outgroup, type_of_scaling=self.variance_correction)/m_scaler(self.bias_c_weight, p, n_outgroup)
                if self.add_variance_correction_to_graph:
                    if self.save_variance_correction:
                        savetxt(self.prefix_for_saving_variance_correction+'variance_correction.txt', changer)
                else:
                    m-=changer
            
        
        return m