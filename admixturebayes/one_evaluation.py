
def one_evaluation(starting_tree, 
                   posterior_function,
                   result_file, 
                   alternative_empirical_matrices=[]):
    likelihood_val, prior_val= posterior_function(starting_tree, verbose=True)
    posterior_val = likelihood_val+prior_val
    likelihood_val2= posterior_function.get_max_likelihood(verbose=True)
    likelihood_val3= posterior_function.get_non_empirical_max_likelihood(starting_tree, verbose=True)
    other_statistics= posterior_function.get_size_diff(starting_tree)
    alt_vals=[]
    for empcovmat in alternative_empirical_matrices:
        alt_vals.append(posterior_function.alternative_emp_cov_likelihood(empcovmat, starting_tree, verbose=True))
    
    with open(result_file, 'w') as f:
        f.write(' '.join(['starting_tree']+list(map(str,[likelihood_val, prior_val, posterior_val])))+'\n')
        f.write(' '.join(['max_likelihood']+list(map(str,[likelihood_val2, '.', '.'])))+'\n')
        f.write(' '.join(['pmax_likelihood']+list(map(str,[likelihood_val3, '.', '.'])))+'\n')
        for name, stat in zip(['median_ratio', 'relative_sign', 'frobenius'], other_statistics):
            f.write(' '.join([name]+list(map(str,[stat, '.', '.'])))+'\n')
        for n, alt_val in enumerate(alt_vals):
            f.write(' '.join(['alternative_'+str(n)]+list(map(str,[alt_val, '.', '.'])))+'\n')
    