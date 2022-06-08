from Rtree_to_covariance_matrix import make_covariance
from reduce_covariance import reduce_covariance
from scipy.stats import wishart, norm
from numpy.linalg import eig, LinAlgError
from numpy import sum, sqrt

def n_mark(cov_mat):
    m=cov_mat.shape[0]
    w,_=eig(cov_mat)
    S=sum(w)
    print(S)
    USS=sum(w*w)
    print(USS)
    return (m-1)*S**2 / ( (m+1)*USS-S**2 )
    
def likelihood_treemix(x, emp_cov, variances, nodes=None, pks={}):
    tree, add= x
    r=emp_cov.shape[0]
    if nodes is None:
        nodes=["s"+str(i) for i in range(1,r+1)]
    par_cov=make_covariance(tree, nodes)
    #removedprin (par_cov, emp_cov, add)
    pks['covariance']=par_cov
    if par_cov is None:
        print('illegal tree')
        return -float('inf')
    try:
        #removedprin emp_cov-add
        #removedprin add
        #removedprin par_cov
        d=sum(norm.logpdf((emp_cov-par_cov-add)/sqrt(variances)))
    except (ValueError, LinAlgError) as e:
        #removedprin "illegal par_cov matrix or to large add"
        #removedprin e
        return -float("inf")
    return d

def likelihood_treemix_from_matrix(matrix, emp_cov, variances,  pks={}):
    pks['covariance']=matrix
    if matrix is None:
        print('illegal tree')
        return -float('inf')
    try:
        #removedprin emp_cov-add
        #removedprin add
        #removedprin par_cov
        d=sum(norm.logpdf((emp_cov-matrix)/sqrt(variances)))
    except (ValueError, LinAlgError) as e:
        #removedprin "illegal par_cov matrix or to large add"
        #removedprin e
        return -float("inf")
    return d

def likelihood(x, emp_cov, b, M=12,nodes=None, collapse_row='', pks={}):
    tree, add= x
    r=emp_cov.shape[0]
    if nodes is None:
        nodes=["s"+str(i) for i in range(1,r+1)]
    par_cov=make_covariance(tree, nodes)
    if par_cov is None:
        print('illegal tree')
        return -float('inf')
    if collapse_row:
        n=len(nodes)-1   #n_outgroup=next((n for n, e in enumerate(nodes_with_outgroup) if e==outgroup))
        par_cov=reduce_covariance(par_cov, n)
    if b is not None:
        par_cov+=b
    #removedprin (par_cov, emp_cov, add)
    pks['covariance']=par_cov
    if par_cov is None:
        print('illegal tree')
        return -float('inf')
    try:
        d=wishart.logpdf(emp_cov, df=M, scale= (par_cov+add)/M)
    except (ValueError, LinAlgError) as e:
        return -float("inf")
    return d

def likelihood_from_matrix(matrix, emp_cov, b, M,  pks={}):
    if b is not None:
        matrix2=matrix+b
    else:
        matrix2=matrix
    pks['covariance']=matrix2
    if matrix2 is None:
        print('illegal tree')
        return -float('inf')
    try:
        d=wishart.logpdf(emp_cov, df=M, scale= matrix2/M)
    except (ValueError, LinAlgError) as e:
        return -float("inf")
    return d
