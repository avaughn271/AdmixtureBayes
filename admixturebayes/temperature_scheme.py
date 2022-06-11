
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