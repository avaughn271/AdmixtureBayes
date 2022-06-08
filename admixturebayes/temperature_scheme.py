from math import exp

class fixed_geometrical(object):
    
    def __init__(self, maxT, no_chains):
        if no_chains==1:
            self.temps=[1.0]
        else:
            self.temps=[maxT**(float(i)/(float(no_chains)-1.0)) for i in range(no_chains)]
    
    def get_temp(self,i):
        return self.temps[i]
    
    def update_temps(self, permut):
        pass
    
class temperature_adapting(object):
    
    def __init__(self, initial_maxT, no_chains, fixed_maxT=False):
        self.count=10
        self.alpha=0.5
        self.average=[0]*no_chains
        self.no_updates=0
        self.maxT=initial_maxT
        self.fixed_maxT=fixed_maxT
        if no_chains==1:
            self.temps=[1.0]
        else:
            self.temps=[initial_maxT**(float(i)/(float(no_chains)-1.0)) for i in range(no_chains)]
     
    def get_temp(self, i):
        return self.temps[i]
    
    def update_temps(self, permut):
        being_crosseds=[0 if all((permut[i]<=j for i in range(0,j+1))) else 1 for j in range(len(permut)-1)]
        self.count+=1
        diffs=[j/i for i,j in zip(self.temps[:-1],self.temps[1:])]
        gamma=0.5/(self.count**self.alpha)
        update_sizes=[exp(gamma*(0-0.234)), exp(gamma*(1-0.234))]
        self.temps=[1.0]
        for d, crossed in zip(diffs, being_crosseds):
            self.temps.append(max(self.temps[-1],self.temps[-1]*d*update_sizes[crossed]))
        new_temps=[]
        if self.fixed_maxT:
            for t in self.temps:
                new_temps.append(t/self.maxT)
        self.no_updates+=1
        self.average=[av+cr for av,cr in zip(self.average, being_crosseds)]
        print('new temperatures:', [round(f,3) for f in self.temps])
        print('average acceptances', [round(float(f)/self.no_updates,3) for f in self.average])
