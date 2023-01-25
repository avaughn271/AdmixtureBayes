from copy import deepcopy

from numpy import insert, delete, ix_, dtype, loadtxt
import subprocess
import os

from numpy import array, mean, zeros, diag, savetxt, nan, isnan, nanmean, identity, diag, log, outer, square
from numpy import sum as npsum
import warnings

from numpy.random import choice
from numpy.linalg import norm as npnorm
from numpy import var as npvar
from scipy.optimize import minimize_scalar

from pathos.multiprocessing import Pool

def gzip(filename, new_filename=None):
    if new_filename is None:
        new_filename=filename+'.gz'
    command=['gzip','-c',filename]
    with open(new_filename, 'w') as f:
        subprocess.call(command, stdout=f)
    return new_filename

def variance_mean_based(sample_of_matrices):
    mean_wishart=  mean(sample_of_matrices, axis=0)
    var_wishart= npvar(sample_of_matrices, axis=0)
    var_rom_mean_wishart=square(mean_wishart)+outer(diag(mean_wishart),diag(mean_wishart))

    def fff(exx):
        return(log(npnorm(var_wishart-var_rom_mean_wishart/exx)**2))
    
    yyy = minimize_scalar(fff, method='bounded', bounds=[0,1000000]).x
    assert yyy > 2.0, "Bootstrap Optimization Failed"
    assert yyy < 999900, "Bootstrap Optimization Failed"
    return yyy
    #maximizes the function, function in the interval [lower_limit,\infty).

def get_partitions(lines, blocksize):
    list_of_lists=[]
    for i in range(0,len(lines)-blocksize, blocksize):
        list_of_lists.append(lines[i:(i+blocksize)])
    if len(lines)%blocksize == 0:
        list_of_lists.append(lines[-blocksize:])
    return list_of_lists

def combine_covs(tuple_covs, indices):
    cov_sum, scale_sum=0.0,0.0
    for i in indices:
        cov_sum+=tuple_covs[i][0]
        scale_sum+=tuple_covs[i][1]
    return cov_sum/scale_sum

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

def make_single_files(filename,blocksize, verbose_level='normal'):
    filenames=[]
    os.mkdir(os.getcwd() + "/temp_adbayes")
    filename_reduced=os.getcwd() + "/temp_adbayes/" +filename.split(os.sep)[-1]+'boot.'
    with open(filename, 'r') as f:
        first_line=f.readline()
        lines=f.readlines()
    n=len(lines)
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
    return filenames

def estimate_degrees_of_freedom_scaled_fast(filename, bootstrap_blocksize=1000,
                                            cores=1,  verbose_level='normal',  **kwargs):
    single_files=make_single_files(filename, blocksize=bootstrap_blocksize, verbose_level=verbose_level)
    assert len(single_files)>1, 'There are ' +str(len(single_files)) + ' bootstrapped SNP blocks and that is not enough. Either add more data or lower the --bootstrap_blocksize'
    if len(single_files)<39:
        warnings.warn('There are only '+str(len(single_files))+' bootstrap blocks. Consider lowering the --bootstrap_blocksize or add more data.', UserWarning)
    single_covs=make_covariances(single_files, cores=cores, return_also_mscale=True, **kwargs)

    p=Pool(cores)
    def t(empty):
        indices=choice(len(single_covs), size=len(single_covs), replace = True)
        return combine_covs(single_covs, indices)
    return variance_mean_based(list(map(t, list(range(100)))))

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
    return mean(array([pi*(1-pi)/(ni-1) for pi,ni in zip(p,n) if ni>1]))

def reduced_covariance_bias_correction(p,n,n_outgroup=0):
    Bs=[]
    for pi,ni in zip(p,n):
        Bs.append(var(pi,ni))
    outgroup_b=Bs.pop(n_outgroup)
    return diag(array(Bs))+outgroup_b

class ScaledEstimator(object):
    def __init__(self, add_variance_correction_to_graph=False, save_variance_correction=True, nodes=None):
        self.add_variance_correction_to_graph=add_variance_correction_to_graph
        self.nodes=nodes
        self.save_variance_correction=save_variance_correction
        
    def __call__(self, xs, ns, extra_info={}):
        if 0 in ns:
            warnings.warn('There were 0s in the allele-totals, inducing nans and slower estimation.', UserWarning)
            ps=nan_divide(xs, ns)
        else:
            ps=xs/ns 
        return self.estimate_from_p(ps, ns=ns, extra_info=extra_info)
        
    def estimate_from_p(self, p, ns=None, extra_info={}):
        p2 = p-p[0,:]
        if npsum(isnan(p2))>0:
            warnings.warn('Nans found in the allele frequency differences matrix => slower execution', UserWarning)
            m=nan_product(p2, p2.T)
        else:
            m=p2.dot(p2.T)/p2.shape[1]
        
        scaling_factor=m_scaler(p)
        extra_info['m_scale']=scaling_factor
        m=m/scaling_factor
        m=reduce_covariance(m)
        b=reduced_covariance_bias_correction(p, ns, 0)/scaling_factor
        if self.save_variance_correction:
            savetxt('variance_correction.txt', b)
        return m

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

def make_estimator( nodes,  reducer, add_variance_correction_to_graph=False, save_variance_correction=True):
    return ScaledEstimator(add_variance_correction_to_graph=add_variance_correction_to_graph, save_variance_correction=save_variance_correction)

def rescale_empirical_covariance(m):
    n=m.shape[0]
    max_expected_trace=n*(n+1)/2-1
    multiplier= max_expected_trace/m.trace()
    
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

    if est_args['save_variance_correction']:
        filename='variance_correction.txt'
        vc=loadtxt(filename)
        vc=reorder_reduced_covariance(vc, names, est_args['nodes'], outgroup=est_args['reducer'])
        savetxt(filename, vc)
    if 'return_also_mscale' in kwargs:
        return cov, extra_info_dic['m_scale']
    return cov

def get_covariance(input, full_nodes=None, reduce_covariance_node=None, estimator_arguments={}):
    kwargs={}
    after_reduce_nodes=deepcopy(full_nodes)
    after_reduce_nodes.remove(reduce_covariance_node)
    kwargs['est']=estimator_arguments

    statistic = input
    statistic = empirical_covariance_wrapper_directly(statistic, **kwargs)
    emp_cov_to_file(statistic, 'covariance_without_reduce_name.txt', after_reduce_nodes)

    statistic=rescale_empirical_covariance(statistic)
    filename='covariance_and_multiplier.txt'
    emp_cov_to_file(statistic[0], filename, after_reduce_nodes)
    with open(filename, 'a') as f:
        f.write('multiplier='+str(statistic[1]))
    return statistic
