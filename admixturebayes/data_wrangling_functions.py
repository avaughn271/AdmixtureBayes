from collections import Counter
from numpy import isfinite

def values_to_numbers(list_of_objects):
    res=[]
    seen={}
    count=0
    for element in list_of_objects:
        if element not in seen:
            count+=1
            seen[element]=count
        res.append(seen[element])
    return res

def count_strings(values):
    dic= Counter(values)
    n=len(values)
    return list(dic.keys()), _normalize(list(dic.values()), n)

def _normalize(values, normalizer):
    return [float(val)/normalizer for val in values]

def count_strings2(values1, values2):
    dic1=Counter(values1)
    dic2=Counter(values2)
    counts1=[]
    counts2=[]
    keys1=list(dic1.keys())
    keys2=list(dic2.keys())
    keys=list(set(keys1+keys2))
    for key in keys:
        counts1.append(dic1.get(key,0))
        counts2.append(dic2.get(key,0))
    return keys,_normalize(counts1, len(values1)), _normalize(counts2, len(values2))

def thin_out_nans(x, index):
    return list(zip(*[(y,i) for y,i in zip(x,index) if isfinite(y)]))