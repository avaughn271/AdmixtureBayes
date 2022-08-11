from copy import deepcopy

from numpy import insert, delete, ix_, dtype, loadtxt
import subprocess
import os

from numpy import array, mean, zeros, diag, savetxt, nan, isnan, nanmean, identity
from numpy import sum as npsum
import warnings

class Estimator(object):
    
    def __init__(self):
        pass
    
    def __call__(self,xs,ns, names=None, extra_info={}):
        '''
        Should return the covariance estimate and possibly set the variable fitted_value
        '''
        pass

def reduce_covariance(covmat):
    reducer=insert(identity(covmat.shape[0]-1), 0, -1, axis=1)
    return reducer.dot(covmat).dot(reducer.T)

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

def var(p,n):
    entries=array([pi*(1-pi)/(ni-1) for pi,ni in zip(p,n) if ni>1])
    return mean(entries)

def reduced_covariance_bias_correction(p,n,n_outgroup=0):
    Bs=[]
    for pi,ni in zip(p,n):
        Bs.append(var(pi,ni))
    outgroup_b=Bs.pop(n_outgroup)
    return diag(array(Bs))+outgroup_b

class ScaledEstimator(Estimator):
    def __init__(self,
                 add_variance_correction_to_graph=False,
                 save_variance_correction=True,
                 nodes=None):
        super(ScaledEstimator, self).__init__()
        self.variance_correction='unbiased'
        self.add_variance_correction_to_graph=add_variance_correction_to_graph
        self.nodes=nodes
        self.save_variance_correction=save_variance_correction
        
    def subtract_ancestral_and_get_outgroup(self,p):
        return p-p[0,:]
        
    def __call__(self, xs, ns, extra_info={}):
        if 0 in ns:
            warnings.warn('There were 0s in the allele-totals, inducing nans and slower estimation.', UserWarning)
            ps=nan_divide(xs, ns)
        else:
            ps=xs/ns 
        return self.estimate_from_p(ps, ns=ns, extra_info=extra_info)
        
    def estimate_from_p(self, p, ns=None, extra_info={}):
        
        p2 = self.subtract_ancestral_and_get_outgroup(p)
        if npsum(isnan(p2))>0:
            warnings.warn('Nans found in the allele frequency differences matrix => slower execution', UserWarning)
            m=nan_product(p2, p2.T)
        else:
            m=p2.dot(p2.T)/p2.shape[1]
        
        scaling_factor=m_scaler( p)
        extra_info['m_scale']=scaling_factor
        m=m/scaling_factor
        m=reduce_covariance(m)
        if self.variance_correction!='None':
            assert ns is not None, 'Variance correction needs a ns-matrix specified'
            b=reduced_covariance_bias_correction(p, ns, 0)/scaling_factor
            if self.save_variance_correction:
                savetxt('variance_correction.txt', b)
        
        return m

def thin_covariance(covmat, nodes_order, specified_nodes):
    ni={node:i for i,node in enumerate(nodes_order)}
    take_out_indices=[ni[s] for s in specified_nodes]
    return covmat[ix_(take_out_indices,take_out_indices)]

def read_freqs(new_filename):
    with open(new_filename, 'r') as f:
        names=f.readline().split()
        pop_sizes=[]
        minors=[]
        for n,r in enumerate(f.readlines()):
            minor_majors=r.split()
            minor_list=[]
            pop_sizes_SNP=[]
            for minor_major in minor_majors:
                minor, major= list(map(float,minor_major.split(',')))
                minor_list.append(minor)
                pop_sizes_SNP.append(major+minor)
            minors.append(minor_list)
            pop_sizes.append(pop_sizes_SNP)
    return names, pop_sizes, minors

def get_xs_and_ns_from_treemix_file(snp_file):
    if snp_file.endswith('.gz'):
        new_filename=unzip(snp_file)
    else:
        new_filename=snp_file
    names, ns, minors= read_freqs(new_filename)
    xs=array(minors, dtype=dtype(float)).T
    ns=array(ns, dtype=dtype(float)).T
    return xs,ns,names

def order_covariance(xnn_tuple, outgroup=''):
    if not outgroup:
        return xnn_tuple
    xs,ns,names=xnn_tuple
    assert outgroup in names, 'The outgroup was not found in the data. Did you spell it correctly?'
    n_outgroup=next((n for n, e in enumerate(names) if e==outgroup))
    xs_o=xs[n_outgroup,:]
    ns_o=ns[n_outgroup,:]
    names_o=names[n_outgroup]
    xs=delete(xs, n_outgroup,0)
    ns=delete(ns, n_outgroup,0)
    names.remove(names_o)
    xs=insert(xs, 0, xs_o, axis=0)
    ns=insert(ns, 0, ns_o, axis=0)
    names=[names_o]+names
    return xs,ns,names

def _get_permutation(actual, target):
    val_to_index={val:key for key, val in enumerate(actual)}
    indices_to_get_target=[val_to_index[val] for val in target]
    return indices_to_get_target

def reorder_reduced_covariance(cov, names, full_nodes, outgroup=''):
    names2=deepcopy(names)
    full_nodes2=deepcopy(full_nodes)
    names2.remove(outgroup)
    full_nodes2.remove(outgroup)
    indices=_get_permutation(names2, full_nodes2)
    return cov[ix_(indices, indices)]

def emp_cov_to_file(m, filename='emp_covimport', nodes=None):
    with open(filename, 'w') as f:
        f.write(' '.join(nodes)+'\n')
        for i, node in enumerate(nodes):
            f.write(node+ ' '+ ' '.join(map(str, m[i]))+'\n')

def unzip(filename, overwrite=False, new_filename=None):
    if new_filename is None:
        new_filename=filename[:-3]
    if (not overwrite) and os.path.exists(new_filename):
        return new_filename
    command=['gunzip','-c',filename]
    with open(new_filename, 'w') as f:
        subprocess.call(command, stdout=f)
    return new_filename

def make_estimator( nodes, 
                   reducer,
                   Simulator_fixed_sxeed,
                   add_variance_correction_to_graph=False,
                   save_variance_correction=True):
    
    est=ScaledEstimator(add_variance_correction_to_graph=add_variance_correction_to_graph,
                 save_variance_correction=save_variance_correction)
    return(est)

def rescale_empirical_covariance(m):
    n=m.shape[0]
    actual_trace=m.trace()
    max_expected_trace=n*(n+1)/2-1
    
    multiplier= max_expected_trace/actual_trace
    
    return m*multiplier, multiplier

def empirical_covariance_wrapper_directly(snp_data_file, **kwargs):
    xnn_tuple=get_xs_and_ns_from_treemix_file(snp_data_file)
    return xnn_to_covariance_wrapper_directly(xnn_tuple, **kwargs)

def xnn_to_covariance_wrapper_directly(xnn_tuple, **kwargs):
    est_args=kwargs['est']
    xnn_tuple=order_covariance(xnn_tuple, outgroup=est_args['reducer'])
    xs,ns,names=xnn_tuple

    est= make_estimator( **est_args)
    extra_info_dic={}
    cov=est(xs,ns, extra_info_dic)
    cov=reorder_reduced_covariance(cov, names, est_args['nodes'], outgroup=est_args['reducer'])
    if ('add_variance_correction_to_graph' in est_args and
        est_args['add_variance_correction_to_graph'] and
        'save_variance_correction' in est_args and
        est_args['save_variance_correction']):
        filename='variance_correction.txt'
        vc=loadtxt(filename)
        vc=reorder_reduced_covariance(vc, names, est_args['nodes'], outgroup=est_args['reducer'])
        savetxt(filename, vc)
    if 'm_scale' in extra_info_dic:
        if 'return_also_mscale' in kwargs and kwargs['return_also_mscale']:
            return cov, extra_info_dic['m_scale']
    return cov

def normaliser_wrapper(covariance, **kwargs):
    return rescale_empirical_covariance(covariance)

dictionary_of_transformations={
    (6,8):empirical_covariance_wrapper_directly,
    (8,9):normaliser_wrapper
    }

dictionary_of_reasonable_names = { 8:'covariance_without_reduce_name', 9:'covariance_and_multiplier'}

def save_stage(value, stage_number, full_nodes, after_reduce_nodes):
    save_word=dictionary_of_reasonable_names[stage_number]
    filename=save_word+'.txt'
    if stage_number==8:
        emp_cov_to_file(value, filename, after_reduce_nodes)
    else:
        emp_cov_to_file(value[0], filename, after_reduce_nodes)
        with open(filename, 'a') as f:
            f.write('multiplier='+str(value[1]))

def get_covariance(input, full_nodes=None,
                   reduce_covariance_node=None,
                   estimator_arguments={}):

    kwargs={}
    after_reduce_nodes=deepcopy(full_nodes)
    after_reduce_nodes.remove(reduce_covariance_node)
    kwargs['est']=estimator_arguments

    stages_to_go_through = [6,8,9]
    statistic = input
    for stage_from, stage_to in zip(stages_to_go_through[:-1], stages_to_go_through[1:]):
        transformer_function=dictionary_of_transformations[(stage_from, stage_to)]
        statistic=transformer_function(statistic, **kwargs)
        save_stage(statistic, stage_to, full_nodes, after_reduce_nodes)
    return statistic
