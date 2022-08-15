from copy import deepcopy

from numpy import insert, delete, ix_, dtype, loadtxt
import subprocess
import os

from numpy import array, mean, zeros, diag, savetxt, nan, isnan, nanmean, identity
from numpy import sum as npsum
import warnings

from numpy.random import choice

from pathos.multiprocessing import Pool

import numpy as np

def gzip(filename, new_filename=None):
    if new_filename is None:
        new_filename=filename+'.gz'
    command=['gzip','-c',filename]
    with open(new_filename, 'w') as f:
        subprocess.call(command, stdout=f)
    return new_filename



#maximizes the function, function in the interval [lower_limit,\infty).
def I_cant_believe_I_have_to_write_this_function_myself(function, lower_limit):
    old_x=lower_limit
    new_x=lower_limit*2
    old_y=function(old_x)
    lower_y=old_y
    new_y=function(new_x)
    sgn=1
    max_step_size_increases=20
    step_size_increases=0
    step_size=lower_limit
    for _ in range(20):
        c=0
        while new_y<old_y and c<100:
            old_x=new_x
            new_x+=sgn*step_size
            old_y=new_y
            if new_x<lower_limit:
                new_x=lower_limit
                new_y=lower_y
                break
            new_y=function(new_x)
            c+=1
            if c>10 and max_step_size_increases>step_size_increases:
                c=0
                step_size*=2
                step_size_increases+=1
        sgn=-sgn
        step_size*=0.5
        old_x=new_x
        new_x=old_x+sgn*step_size
        old_y=new_y
        new_y=function(new_x)
    return new_x

def variance_mean_based(sample_of_matrices, divisor=None, verbose_level='normal'):
    mean_wishart=  np.mean(sample_of_matrices, axis=0)
    var_wishart= np.var(sample_of_matrices, axis=0)
    r=mean_wishart.shape[0]
    var_rom_mean_wishart=np.square(mean_wishart)+np.outer(np.diag(mean_wishart),np.diag(mean_wishart))
    def penalty_function(df_l):
        df=df_l
        val=np.linalg.norm(var_wishart-var_rom_mean_wishart/df)**2
        return np.log(val)
    
    rval=I_cant_believe_I_have_to_write_this_function_myself(penalty_function, r)
    return rval

def get_partitions(lines, blocksize):
    list_of_lists=[]
    for i in range(0,len(lines)-blocksize, blocksize):
        list_of_lists.append(lines[i:(i+blocksize)])
    if len(lines)%blocksize == 0:
        list_of_lists.append(lines[-blocksize:])
    return list_of_lists

def bootstrap_indices(k):
    bootstrap_inds=choice(k, size=k, replace = True)
    return bootstrap_inds

def combine_covs(tuple_covs, indices):
    cov_sum, scale_sum=0.0,0.0
    for i in indices:
        cov_sum+=tuple_covs[i][0]
        scale_sum+=tuple_covs[i][1]
    return cov_sum/scale_sum

def bootsrap_combine_covs(covs, cores, bootstrap_samples):
    p=Pool(cores)
    def t(empty):
        indices=bootstrap_indices(len(covs))
        return combine_covs(covs, indices)
    result_covs=list(map(t, list(range(bootstrap_samples))))
    return result_covs

def remove_files(filenames):
    for fil in filenames:
        os.remove(fil)
        os.remove(fil[:-3])

def make_covariances(filenames, cores, **kwargs):
    covs=[]
    p=Pool(cores)
    def t(filename):
        return empirical_covariance_wrapper_directly(filename, **kwargs)
    covs=list(map(t, filenames))
    try:
        remove_files(filenames)
    except OSError as e:
        warnings.warn('Erasing the files did not succeed',UserWarning)
    return covs

def make_single_files(filename,blocksize, no_blocks, verbose_level='normal'):
    filenames=[]
    os.mkdir(os.getcwd() + "/temp_adbayes")
    filename_reduced=os.getcwd() + "/temp_adbayes/" +filename.split(os.sep)[-1]+'boot.'
    with open(filename, 'r') as f:
        first_line=f.readline()
        lines=f.readlines()
    n=len(lines)
    if no_blocks is not None:
        blocksize=n/no_blocks
    line_sets=get_partitions(lines, blocksize)
    if verbose_level!='silent':
        print('total number of SNPs: '+ str(n))
    for i,lins in enumerate(line_sets):
        new_filename= filename_reduced+str(i)
        with open(new_filename, 'w') as g:
            g.write(first_line)
            g.writelines(lins)
        gzipped_filename=gzip(new_filename)
        filenames.append(gzipped_filename)
    return filenames, first_line.split()

def estimate_degrees_of_freedom_scaled_fast(filename,
                                            bootstrap_blocksize=1000,
                                            no_blocks=None,
                                            no_bootstrap_samples=10,
                                            cores=1,
                                            verbose_level='normal',
                                            **kwargs):
    single_files, nodes=make_single_files(filename, blocksize=bootstrap_blocksize, no_blocks=no_blocks, verbose_level=verbose_level)
    assert len(single_files)>1, 'There are ' +str(len(single_files)) + ' bootstrapped SNP blocks and that is not enough. Either add more data or lower the --bootstrap_blocksize'
    if len(single_files)<39:
        warnings.warn('There are only '+str(len(single_files))+' bootstrap blocks. Consider lowering the --bootstrap_blocksize or add more data.', UserWarning)
    single_covs=make_covariances(single_files, cores=cores, return_also_mscale=True, **kwargs)
    covs=bootsrap_combine_covs(single_covs, cores=cores, bootstrap_samples=no_bootstrap_samples)
    res=variance_mean_based(covs, verbose_level=verbose_level)
    return res, covs

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
