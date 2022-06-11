
def get_number_of_admixes(tree):
    tot_length=sum(map(len, tree))
    no_admixs=(tot_length-len(tree)*3)/4
    return no_admixs