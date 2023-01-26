import torch

class SPH_CAL():
    def __init__(self,max_l=10,device=torch.device("cpu"),dtype=torch.float32):
        '''
         max_l: maximum L for angular momentum
        device: cpu/gpu
         dtype:  torch.float32/torch.float64

        '''
        #  form [0,max_L]
        if max_l<1: raise ValueError("The angular momentum must be greater than or equal to 1. Or the angular momentum is lack of angular information, the calculation of the sph is meanless.")
        self.max_l=max_l+1
        self.pt=torch.empty((self.max_l,self.max_l),device=device,dtype=torch.int32)
        self.yr=torch.empty((self.max_l,self.max_l),device=device,dtype=torch.int32)
        self.yr_rev=torch.empty((self.max_l,self.max_l),device=device,dtype=torch.int32)
        num_lm=((self.max_l+1)*self.max_l/2-2
        self.coeff_a=torch.empty(num_lm,device=device,dtype=dtype)
        self.coeff_b=torch.empty(num_lm,device=device,dtype=dtype)
        tmp=torch.arange(self.max_l).to(device).to(dtype)
        ls=tmp*tmp
        lms=(tmp-1)*(tmp-1)
        for m in range(self.max_l):
            self.pt[m,m:self.max_l]=m+tmp[m:self.max_l]*(tmp[m:self.max_l]+1)/2
            self.yr[m,m:self.max_l]=m+tmp[m:self.max_l]*(tmp[m:self.max_l]+1)
            self.yr_rev[m,m:self.max_l]=-m+tmp[m:self.max_l]*(tmp[m:self.max_l]+1)
            if m<self.maxl_l-2:
                self.coeff_a[self.pt[m,max(2,m):self.max_l]]=np.sqrt((4.0*ls[max(2,m):self.max_l]-1.0)/(ls[max(2,m):self.max_l]-lms[max(2,m):self.max_l]))
                self.coeff_b[self.pt[m,max(2,m):self.max_l]]=-np.sqrt((lmls[max(2,m):self.max_l]-m*m)/(4.0*lmls[max(2,m):self.max_l]-1.0))

    def __call__(self,cart)
        return self.computeY(cart)

    def computeY(self,cart):
        '''
        cart: Cartesian coordinates with the dimension (3,n) n is the number of total atoms.
        '''
        r=torch.linalg.vector_norm(cart,dim=0)
        cos_ceta=cart[2,:]/r
        sin_ceta=torch.sqrt(1.0-cos_ceta*cos_ceta)
        ALP=self.computeP(self,cos_ceta,sin_ceta)
        sqrt2_div2=np.sqrt(2)/2.0
        sph=torch.empty((self.max_l*self.max_l,cart.shape[1]),device=cart.device,dtype=cart.dtype)
        sph[self.yr[0,0:max_l]]=ALP[self.pt[0,0:max_l]]*sqrt2_div2

        r1=r*sin_ceta
        c1=1.0
        c2=cart[0,:]/r1
        s1=0.0
        s2=-cart[1,:]/r1
        dc2=2.0*c2
        for m in range(1,self.max_l):
            s=dc2*s1-s2
            c=dc2*c1-c2
            s2=s1
            s1=s
            c2=c1
            c1=c
            sph[self.yr_rev[m,m:self.max_l]]=ALP[self.pt[m,m:self.max_l]]*s[None,:]
            sph[self.yr[m,m:self.max_l]]=ALP[self.pt[m,m:self.max_l]]*c[None,:]

    def computeP(self,cos_ceta,sin_ceta)
        temp=np.sqrt(0.5/np.pi)
        ALP=torch.empty((self.max_l*(self.max_l+1)/2,cos_ceta.shape[0])device=cos_ceta.device,dtype=cos_ceta.dtype)
        ALP[0]=temp
        ALP[1]=cos_ceta*np.sqrt(3)*temp
        ALP[2]==-np.sqrt(3/2)*sin_ceta*temp
        for l in range(2,self.max_l):
            for m in range(l-1):
                ALP[self.pt[m,l]]=self.coeff_a[self.pt[m,l]]*(cos_ceta*ALP[self.pt[m,l-1]]+self.coeff_b[pt[m,l]]*ALP[self.pt[m,l-2]])
            ALP[self.pt[l-1,l]]=cos_ceta*np.sqrt(2*l+1)*ALP[self.pt[l-1,l-1]]
            ALP[self.pt[l,l]]=-np.sqrt(1.0+1.0/2.0/l)*sin_ceta*ALP[self.pt[l-1,l-1]]
        return ALP

