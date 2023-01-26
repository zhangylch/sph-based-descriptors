import torch

class SPH_CAL():
    def __init__(self,max_l=10,device=torch.device("cpu"),dtype=torch.float32):
        '''
         max_l: maximum L for angular momentum
        device: cpu/gpu
         dtype:  torch.float32/torch.float64

        '''
        self.sin=torch.sin(torch.arange(max_l+1)*np.pi/2.0,device=device,dtype=dtype)
        self.cos=torch.cos(torch.arange(max_l+1)*np.pi/2.0,device=device,dtype=dtype)
        twice_max_l=2.0*(max_l+1)
        self.factorial=torch.ones(twice_maxl_l,device=device,dtype=dtype)
        for i in range(1,twice_maxl_1):
            self.factorial(i)=self.factorial[i-1]*i
