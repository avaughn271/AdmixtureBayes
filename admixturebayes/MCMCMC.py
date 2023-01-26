import pandas as pd
from MCMC import basic_chain_pool
from itertools import chain

def MCMCMC(coolerMultiple, starting_trees,    posterior_function, summaries, temperature_scheme, 
           iteration_scheme,  proposal_scheme, n_arg, verboseee,
           multiplier= None, result_file=None):
    df_result=None
    xs = starting_trees

    pool = basic_chain_pool(summaries, posterior_function, proposal_scheme)
    posteriors = [posterior_function(x) for x in xs]

    proposal_updates=[proposal.get_exportable_state() for proposal in proposal_scheme]

    cum_iterations=0
    print(n_arg)
    for no_iterations in iteration_scheme:
        if cum_iterations % 1000 == 0 and verboseee != "silent":
            print("Currently on iteration " +  str(cum_iterations) + " out of " + str(n_arg * 50))
        iteration_object=_pack_everything(xs, posteriors, temperature_scheme, no_iterations, cum_iterations, proposal_updates, multiplier)
        new_state = pool.order_calculation(iteration_object)

        xs, posteriors, df_add, proposal_updates = _unpack_everything(new_state, summaries)
        df_result=_update_results(df_result, df_add)
        if result_file is not None:
            if cum_iterations==0:
                start_data_frame(df_result, result_file)
            elif df_result.shape[0]>1000:
                add_to_data_frame(df_result, result_file)
                df_result=df_result[0:0]
        cum_iterations+=no_iterations
        temperature_scheme[0] = temperature_scheme[0] * coolerMultiple
    for chain in pool.group:
            chain.process.terminate()

def start_data_frame(df, result_file):
    df=df.loc[df.layer==0,:]
    df.to_csv(result_file, header=True)

def add_to_data_frame(df_add, result_file):
    df_add=df_add.loc[df_add.layer==0,:]
    with open(result_file, 'a') as f:
        df_add.to_csv(f, header=False)

def _update_results(df_result, df_add):
    if df_result is None:
        df_result = df_add
    else:
        df_result = pd.concat([df_result, df_add])
    return df_result

def _pack_everything(xs, posteriors, temperature_scheme,no_iterations,cum_iterations, proposal_updates=None, multiplier=None):
    return ([x, posterior,no_iterations, 40, 
    cum_iterations, temperature_scheme[i], proposal_update, multiplier] for i,(x, 
    posterior,proposal_update) in enumerate(zip(xs,posteriors,proposal_updates)))

def _unpack_everything(new_state, summaries):
    xs,posteriors, summs, proposal_updates = list(zip(*new_state))
    list_of_smaller_data_frames=[]
    for summ_data, i in zip(summs, list(range(1))):
        iter_chain=chain((('iteration', summ_data[0]),),
                         ((summ_object.name,summ_col) for summ_object,summ_col in zip(summaries,summ_data[1:])))
        df=pd.DataFrame.from_dict(dict(iter_chain))
        df['layer']=i
        list_of_smaller_data_frames.append(df)
    df=pd.concat(list_of_smaller_data_frames)
    return list(xs), list(posteriors), df, list(proposal_updates)
