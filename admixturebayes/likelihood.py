from Rtree_to_covariance_matrix import make_covariance
from reduce_covariance import reduce_covariance
from scipy.stats import wishart
from numpy.linalg import LinAlgError

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
        n=len(nodes)-1
        par_cov=reduce_covariance(par_cov, n)
    if b is not None:
        par_cov+=b
    pks['covariance']=par_cov
    if par_cov is None:
        print('illegal tree')
        return -float('inf')
    try:
        d=wishart.logpdf(emp_cov, df=M, scale= (par_cov+add)/M)
    except (ValueError, LinAlgError) as e:
        return -float("inf")
    return d
