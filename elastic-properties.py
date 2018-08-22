#!/usr/bin/env python
# coding: utf-8

# ## Elastic Properties
# ### Code for extracting elastic tensor and calculating mechanical properties from VASP OUTCAR
# 
# Equations can be found at https://www.materialsproject.org/wiki/index.php/Elasticity_calculations

# In[86]:

import numpy as np
import sys


# In[87]:

def get_elastic_tensor(filename):
    ''' Reads the elastic tensor from the OUTCAR. 
    Args:
        filename : the name of the vasp OUTCAR
    Returns:
        elastic_tensor : 6x6 tensor of the elastic moduli
    '''
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    copy = False
    elastic_tensor = []
    for line in lines:
        inp = line.split()
        if inp == []:
            continue 
        if len(inp) < 4 or len(inp) > 7:
            continue
        if len(inp) == 4 and inp[0] == 'TOTAL':
            copy = True
        if copy:
            if len(inp) == 7 and len(inp[0]) == 2:
                elastic_tensor.append(inp[1:])
    return np.asarray(elastic_tensor).astype(np.float)


# ### Elastic tensor $C_{ij}$

# In[88]:


elastic_tensor = get_elastic_tensor('OUTCAR')


# ### Divide by 10 to convert kBar to GPa

# In[89]:

Cij = elastic_tensor/10
print "\nCij="
print Cij

# ### Compliance tensor $s_{ij}$ $(GPa^{-1})$
# $s_{ij} = C_{ij}^{-1}$

# In[90]:

Sij = np.linalg.inv(Cij)
print "\nSij="
print Sij


# ### Voigt bulk modulus $K_v$ $(GPa)$

# In[91]:

Kv = ((Cij[0,0] + Cij[1,1] + Cij[2,2]) + 2 * (Cij[0,1] + Cij[1,2] + Cij[2,0])) / 9


# ### Reuss bulk modulus $K_R$ $(GPa)$

# In[92]:

Kr = 1/((Sij[0,0] + Sij[1,1] + Sij[2,2]) + 2 * (Sij[0,1] + Sij[1,2] + Sij[2,0])) 


# ### Voigt shear modulus $G_v$ $(GPa)$

# In[93]:

Gv = ((Cij[0,0] + Cij[1,1] + Cij[2,2] - Cij[0,1] - Cij[0,2] - Cij[1,2]) / 15) + ((Cij[3,3] + Cij[4,4] + Cij[5,5]) / 5)


# ### Reuss shear modulus $G_v$ $(GPa)$

# In[94]:

Gr = 15 / (4 * (Sij[0,0] + Sij[1,1] + Sij[2,2]) - 4 * (Sij[0,1] + Sij[1,2] + Sij[2,0]) + 3 * (Sij[3,3] + Sij[4,4] + Sij[5,5]))


# ### Voigt-Reuss-Hill bulk modulus $K_{VRH}$ $(GPa)$

# In[95]:

Kvrh = (Kv + Kr)/2
print "Voigt-Reuss-Hill bulk modulus (GPa): ", Kvrh

print "Voigt bulk modulus (GPa): ", Kv
print "Reuss bulk modulus (GPa): ", Kr


# ### Voigt-Reuss-Hill shear modulus $G_{VRH}$ $(GPa)$
# $G_{VRH} = (G_R + G_v)/2$

# In[96]:

Gvrh = (Gv + Gr)/2
print "Voigt-Reuss-Hill shear modulus (GPa): ", Gvrh
print "Voigt shear modulus (GPa): ", Gv
print "Reuss shear modulus (GPa): ", Gr


# ### Isotropic Poisson ratio $\mu$
# $\mu = (3K_{VRH} - 2G_{VRH})/(6K_{VRH} + 2G_{VRH})$

# In[83]:

mu = (3 * Kvrh - 2 * Gvrh) / (6 * Kvrh + 2 * Gvrh )
print "Isotropic Poisson ratio: ", mu

ye = (9* Kvrh * Gvrh) / (3 * Kvrh + Gvrh)
print " young's modulus (GPa)" , ye

yer = (9* Kr * Gr) / (3 * Kr + Gr)
print " test ", yer
# In[85]:

Cij


# In[ ]:




