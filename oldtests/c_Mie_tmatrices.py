#!/usr/bin/env python
# coding: utf-8

# In[1]:


from qpms import *


# In[6]:


R = 40e-9
ω_p = 9*eV/ℏ #9*eV/ℏ
ε_inf = 4.6
γ_p = 0.1*eV/ℏ
ε_b = 2.13
lMax = 3

ω = 1.5*eV/ℏ


# In[7]:


ε_m = ε_drude(ε_inf, ω_p, γ_p, ω)


# In[9]:


k_i = cmath.sqrt(ε_m)*ω/c
k_e = cmath.sqrt(ε_b)*ω/c
RH, RV, TH, TV = mie_coefficients(a=R, nmax=lMax, k_i=k_i, k_e=k_e, J_ext=1, J_scat=3)


# In[11]:


spec = BaseSpec(lMax=lMax)
cT = CTMatrix.spherical(spec, R, k_i, k_e, 1, 1)


# In[16]:


print(np.diag(cT[...]))


# In[18]:


print(RV)
print(RH)


# In[ ]:




