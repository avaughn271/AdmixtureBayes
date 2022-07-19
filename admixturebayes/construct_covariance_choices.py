from tree_statistics import identifier_to_tree_clean, unique_identifier_and_branch_lengths, generate_predefined_list_string
from tree_to_data import (file_to_emp_cov,
                         emp_cov_to_file,
                          get_xs_and_ns_from_treemix_file, order_covariance, reorder_reduced_covariance)
from scipy.stats import wishart
from copy import deepcopy
from reduce_covariance import rescale_empirical_covariance
from numpy import loadtxt, savetxt
from construct_estimator_choices import make_estimator

def empirical_covariance_wrapper_directly(snp_data_file, **kwargs):
    xnn_tuple=get_xs_and_ns_from_treemix_file(snp_data_file, kwargs['locus_filter'])
    return xnn_to_covariance_wrapper_directly(xnn_tuple, **kwargs)

def xnn_to_covariance_wrapper_directly(xnn_tuple, **kwargs):
    est_args=kwargs['est']
    xnn_tuple=order_covariance(xnn_tuple, outgroup=est_args['reducer'])
    xs,ns,names=xnn_tuple

    est= make_estimator(reduce_method='outgroup',
                   reduce_also=True,
                   ns=ns,**est_args)
    extra_info_dic={}
    cov=est(xs,ns, extra_info_dic)
    cov=reorder_reduced_covariance(cov, names, est_args['nodes'], outgroup=est_args['reducer'])
    if ('add_variance_correction_to_graph' in est_args and
        est_args['add_variance_correction_to_graph'] and
        'save_variance_correction' in est_args and
        est_args['save_variance_correction']):
        filename=est_args['prefix']+'variance_correction.txt'
        vc=loadtxt(filename)
        vc=reorder_reduced_covariance(vc, names, est_args['nodes'], outgroup=est_args['reducer'])
        savetxt(filename, vc)
    if 'm_scale' in extra_info_dic:
        if 'import' in kwargs:
            with open(kwargs['import'], 'w') as f:
                txt=str(extra_info_dic['m_scale'])
                assert '\0' not in txt, 'binary content in m_scale file'
                f.write(txt)
        if 'return_also_mscale' in kwargs and kwargs['return_also_mscale']:
            return cov, extra_info_dic['m_scale']
    return cov

def add_wishart_noise(matrix, df):
    r=matrix.shape[0]
    m=wishart.rvs(df=df, scale=matrix/df)
    return m

def normaliser_wrapper(covariance, **kwargs):
    return rescale_empirical_covariance(covariance, normalizer=kwargs['scale_goal'])

dictionary_of_transformations={
    (6,8):empirical_covariance_wrapper_directly,
    (8,9):normaliser_wrapper
    }

dictionary_of_reasonable_names={
    1:'number_of_leaves',
    2:'leaves_admixtures',
    3:'true_tree',
    4:'true_tree_with_outgroup',
    5:'scaled_true_tree',
    6:'SNP_data',
    7:'covariance',
    8:'covariance_without_reduce_name',
    9:'covariance_and_multiplier'}

def write_one_line_to_file(filename, value):
    with open(filename,'w') as f:
        f.write(value)

def write_two_lines_to_file(filename, value1, value2):
    with open(filename, 'w') as f:
        f.write(value1+'\n'+value2)

def save_stage(value, stage_number, prefix, full_nodes, before_added_outgroup_nodes, after_reduce_nodes, filename=None):
    if filename is None:
        save_word=dictionary_of_reasonable_names[stage_number]
        filename=prefix+save_word+'.txt'
    if stage_number==1:
        write_one_line_to_file(filename, str(value))
    elif stage_number==2:
        write_one_line_to_file(filename, str(value))
    elif stage_number==3:
        write_two_lines_to_file(filename, ' '.join(before_added_outgroup_nodes), unique_identifier_and_branch_lengths(value, before_added_outgroup_nodes))
    elif stage_number==4:
        write_two_lines_to_file(filename, ' '.join(full_nodes), unique_identifier_and_branch_lengths(value, full_nodes))
    elif stage_number==5:
        write_two_lines_to_file(filename, ' '.join(full_nodes), unique_identifier_and_branch_lengths(value, full_nodes))
    elif stage_number==6:
        pass
    elif stage_number==7:
        emp_cov_to_file(value, filename, full_nodes)
    elif stage_number==8:
        emp_cov_to_file(value, filename, after_reduce_nodes)
    elif stage_number in [21,22,23]:
        print('Stage not saved')
    else:
        emp_cov_to_file(value[0], filename, after_reduce_nodes)
        with open(filename, 'a') as f:
            f.write('multiplier='+str(value[1]))


def get_covariance(stages_to_go_through, input, full_nodes=None,
                   skewed_admixture_prior_sim=False,
                   p=0.5,
                   outgroup_name=None,
                   add_wishart_noise_to_covariance=False,
                   df_of_wishart_noise_to_covariance=1000,
                   reduce_covariance_node=None,
                   sample_per_pop=50, nreps=2,
                   theta=0.4, sites=500000, recomb_rate=1,
                   ms_file=None,
                   treemix_file=None,
                   blocksize_empirical_covariance=100,
                   scale_tree_factor=0.05,
                   save_stages=list(range(1,6))+list(range(7,10)),
                   prefix='tmp',
                   t_adjust_tree=False,
                   final_pop_size=100.0,
                   via_treemix=True,
                   sadmix=False,
                   scale_goal='min',
                   favorable_init_brownian=False,
                   unbounded_brownian=False,
                   filter_on_outgroup=False,
                   locus_filter=None,
                   estimator_arguments={},
                   verbose_level='normal'):

    if ms_file is None:
        ms_file=prefix+'ms.txt'

    if treemix_file is None:
        treemix_file=prefix+'treemix'

    kwargs={}
    kwargs['skewed_admixture_prior_sim']=skewed_admixture_prior_sim
    kwargs['p']=p
    kwargs['outgroup_name']=outgroup_name
    kwargs['add_wishart_noise_to_covariance']=add_wishart_noise_to_covariance
    kwargs['df_of_wishart_noise_to_covariance']=df_of_wishart_noise_to_covariance
    kwargs['full_nodes']=full_nodes
    kwargs['sadmix']=sadmix
    before_added_outgroup_nodes=deepcopy(full_nodes)
    after_reduce_nodes=deepcopy(full_nodes)
    if outgroup_name is not None and outgroup_name in before_added_outgroup_nodes:
        before_added_outgroup_nodes.remove(outgroup_name)
    if reduce_covariance_node is not None and reduce_covariance_node in after_reduce_nodes:
        after_reduce_nodes.remove(reduce_covariance_node)
    kwargs['reduce_covariance_node']=reduce_covariance_node
    kwargs['after_reduce_nodes']=after_reduce_nodes
    kwargs['before_added_outgroup_nodes']=before_added_outgroup_nodes
    kwargs['sample_per_pop']=sample_per_pop
    kwargs['nreps']=nreps
    kwargs['theta']=theta
    kwargs['sites']=sites
    kwargs['recomb_rate']=recomb_rate
    kwargs['ms_file']=ms_file
    kwargs['treemix_file']=treemix_file
    kwargs['blocksize_empirical_covariance']=blocksize_empirical_covariance
    kwargs['scale_tree_factor']=scale_tree_factor
    kwargs['pks']={}
    kwargs['time_adjust']=t_adjust_tree
    kwargs['final_pop_size']=final_pop_size
    kwargs['via_treemix']=via_treemix
    kwargs['add_file']=prefix+'true_add.txt'
    kwargs['import']=prefix+'m_scale.txt'
    kwargs['scale_goal']=scale_goal
    kwargs['favorable_init_brownian']=favorable_init_brownian
    kwargs['unbounded_brownian']=unbounded_brownian
    kwargs['filter_on_outgroup']=filter_on_outgroup
    kwargs['est']=estimator_arguments
    kwargs['locus_filter']=locus_filter

    #makes a necessary transformation of the input(if the input is a filename or something).
    statistic=read_input(stages_to_go_through[0], input, full_nodes, before_added_outgroup_nodes, after_reduce_nodes)

    if stages_to_go_through[0] in save_stages:
        save_stage(statistic, stages_to_go_through[0], prefix, full_nodes, before_added_outgroup_nodes, after_reduce_nodes)
    for stage_from, stage_to in zip(stages_to_go_through[:-1], stages_to_go_through[1:]):
        print(stage_from, stage_to)
        transformer_function=dictionary_of_transformations[(stage_from, stage_to)]
        statistic=transformer_function(statistic, **kwargs)
        if stage_to in save_stages:
            save_stage(statistic, stage_to, prefix, full_nodes, before_added_outgroup_nodes, after_reduce_nodes)

    return statistic

def write_output(stage, output):
    pass

def read_input(stage, input, full_nodes, before_added_outgroup_nodes, after_reduce_nodes):
    if stage==1:
        return int(input)
    if stage==2:
        return list(map(int, input[1:-1].split(',')))
    if stage==3:
        return read_tree(input, before_added_outgroup_nodes)
    if stage==4:
        return read_tree(input, full_nodes)
    if stage==5:
        return read_tree(input, full_nodes)
    if stage==6:
        return input
    if stage==7:
        return read_covariance_matrix(input, full_nodes)
    if stage==8:
        return read_covariance_matrix(input, after_reduce_nodes)
    if stage==9 : #assuming that it comes in the pair (matrix or matrix filename, multiplier)
        return (read_covariance_matrix(input, after_reduce_nodes), read_multiplier(input))
    assert False, 'The beginning state '+str(stage)+' is unknown.'

def read_multiplier(input):
    with open(input, 'r') as f:
        last_line=f.readlines()[-1]
        return float(last_line.split("=")[1])

def read_covariance_matrix(input, nodes):
    if isinstance(input, str):
        return file_to_emp_cov(input, nodes=nodes)
    else:
        return input

def read_tree(input, nodes):
    if isinstance(input, str):
        if not ';' in input:
            input=read_one_line_skip(filename=input)
            return identifier_to_tree_clean(input, leaves=generate_predefined_list_string(deepcopy(nodes)))
        else:
            return identifier_to_tree_clean(input, leaves=generate_predefined_list_string(deepcopy(nodes)))
    else:
        return input

def read_one_line_skip(filename):
    with open(filename, 'r') as f:
        lines=f.readlines()
        if len(lines[-1])>3:
            return lines[-1].rstrip()
        else:
            return lines[-2].rstrip()

def read_one_line(filename):
    with open(filename, 'r') as f:
        res=f.readline().rstrip()
    return res
