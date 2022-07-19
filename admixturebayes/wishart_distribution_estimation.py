from numpy import var, savetxt, loadtxt
from numpy.random import choice
from tree_to_data import unzip, gzip
from covariance_estimator import initor
import warnings

from construct_covariance_choices import empirical_covariance_wrapper_directly
from pathos.multiprocessing import Pool
from df_estimators import variance_mean_based
import os
    
def estimate(sample_of_matrices):
    return var(sample_of_matrices, axis=0)

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

def make_single_files(filename,blocksize, no_blocks, prefix='', verbose_level='normal'):
    assert (blocksize is not None) or (no_blocks is not None), 'Has to specify either block size or number of blocks'
    filenames=[]
    if filename.endswith('.gz'):
        filename=unzip(filename)
    filename_reduced=prefix+filename.split(os.sep)[-1]+'boot.'
    with open(filename, 'r') as f:
        first_line=f.readline()
        lines=f.readlines()
    n=len(lines)
    if no_blocks is not None:
        blocksize=n/no_blocks
    line_sets=get_partitions(lines, blocksize)
    if verbose_level!='silent':
        print('total number of SNPs: '+ str(n))
        #print('total blocksize: ' + str(blocksize))
    for i,lins in enumerate(line_sets):
        new_filename= filename_reduced+str(i)
        with open(new_filename, 'w') as g:
            g.write(first_line)
            g.writelines(lins)
        gzipped_filename=gzip(new_filename, overwrite=True)
        filenames.append(gzipped_filename)
    return filenames, first_line.split()

def estimate_degrees_of_freedom_scaled_fast(filename, 
                                            bootstrap_blocksize=1000,
                                            no_blocks=None,
                                            no_bootstrap_samples=10,
                                            summarization=['mle_opt','var_opt','opt'],
                                            cores=1,
                                            save_covs='',
                                            prefix='',
                                            load_bootstrapped_covariances=[],
                                            verbose_level='normal',
                                            **kwargs):
    if not load_bootstrapped_covariances:
        single_files, nodes=make_single_files(filename, blocksize=bootstrap_blocksize, no_blocks=no_blocks, prefix=prefix, verbose_level=verbose_level)
        assert len(single_files)>1, 'There are ' +str(len(single_files)) + ' bootstrapped SNP blocks and that is not enough. Either add more data or lower the --bootstrap_blocksize'
        if len(single_files)<39:
            warnings.warn('There are only '+str(len(single_files))+' bootstrap blocks. Consider lowering the --bootstrap_blocksize or add more data.', UserWarning)
        single_covs=make_covariances(single_files, cores=cores, return_also_mscale=True, **kwargs)
        covs=bootsrap_combine_covs(single_covs, cores=cores, bootstrap_samples=no_bootstrap_samples)
        if save_covs:
            for i,cov in enumerate(covs):
                filn=prefix+save_covs+str(i)+'import'
                savetxt(filn, cov)
    else:
        covs=[loadtxt(bcov) for bcov in load_bootstrapped_covariances]
    summarization=initor(summarization)
    res=variance_mean_based(covs, verbose_level=verbose_level)
    return res, covs
