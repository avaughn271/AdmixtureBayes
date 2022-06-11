from math import log, exp
from scipy.stats import uniform


def rvs(fro=0.0, to=1.0):
    guess=uniform.rvs(fro, to-fro)
    while uniform.rvs() > (to-(guess-fro))/(to-fro):
        guess=uniform.rvs(fro, to-fro)
    return guess

def logpdf(x, fro=0.0, to=1.0):
    log_max_height=log((to-fro)*2)
    log_frac=log( (to-(x-fro))/(to-fro) )
    return log_max_height+log_frac
