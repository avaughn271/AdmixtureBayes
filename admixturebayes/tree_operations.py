#this file contains functions that work on trees, mostly to extract information.
def get_tree_flatter_list():
    return [["r","s1",0.3,[">","s3",None,32123], 0.1], 
                       ["s2s3","s2",0.2],
                       ["s2s3","s3",0.15, ["<","s1",0.44,32123],0.05],
                       ["r","s2s3",0.2]]

def get_number_of_admixes(tree):
    tot_length=sum(map(len, tree))
    no_admixs=(tot_length-len(tree)*3)/4
    return no_admixs

def get_number_of_admixes_on_branch(branch):
    return (len(branch)-3)/2

def get_index_of_nonroot_branches(tree):
    res=[-1,-1]
    count=0
    for n,branch in enumerate(tree):
        if branch[0]!="r":
            res[count]=n
            count+=1
            if count==2:
                break
    return [n for n in range(len(tree)) if n not in res]

def get_all_branch_lengths(tree):
    branch_lengths={}
    for branch in tree:
        times=branch[2::2]
        branch_lengths[branch[0]+branch[1]]=sum(times)
    return branch_lengths

def tree_prune(tree, remove_code):
    pruned=[]
    not_pruned=[]
    orphans=[]
    for branch in tree:
        if branch[1] in remove_code:
            pruned.append(branch)
        else:
            if remove_code in branch[0]:
                branch[0]=branch[0].replace(remove_code, '')
            if remove_code in branch[1]:
                branch[1]=branch[1].replace(remove_code, '')
            
            for admixture_event in branch[3::2]:
                if remove_code in admixture_event[0] and remove_code!=admixture_event[0]:
                    admixture_event[0]=admixture_event[0].replace(remove_code, '')
            if branch[0]==branch[1]:
                orphans.append(branch)
            else:
                not_pruned.append(branch)
    return pruned, not_pruned, orphans

def make_flat_list_no_admix(size=2):
    step_size=1.0/size
    res=[["s1s2","s1",step_size], ["s1s2","s2",step_size]]
    rooter="s1s2"
    for i in range(3,size+1):
        old_rooter=rooter
        si="s"+str(i)
        rooter+=si
        res.append([rooter,si, step_size*(i-1)])
        res.append([rooter,old_rooter,step_size])
    res[-2][0]="r"
    res[-1][0]="r"
    return res

def simple_branch_print(branches):
    for branch in branches:
        print(branch)

def pretty_print_branches(branches):
    res=''
    for branch in branches:
        if len(branch)==3:
            res+= str(branch)+'\n'
        else:
            res+= '['+str(branch[0])+', '+str(branch[1])+', '+str(branch[2])+','
            for admixture_event, b_length in zip(branch[3::2],branch[4::2]):
                if admixture_event[2] is None:
                    res+= '\n'+ '  '+str(admixture_event[:2]+admixture_event[3:])+', '+str(b_length)
                else:
                    res+= '\n'+ '  '+str(admixture_event)+', '+str(b_length)
            res+=']'+'\n'
    return res.rstrip()
            

def get_number_of_leaves(tree):
    pass

def illegal_branch_length(tree):
    for branch in tree:
        times=branch[2::2]
        for t in times:
            if t<0:
                return True
    return False
