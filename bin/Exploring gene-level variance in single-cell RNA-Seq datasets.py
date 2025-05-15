#!/usr/bin/env python
# coding: utf-8

# In[100]:


import pandas as pd
import numpy as np


# In[101]:


sc_df = pd.read_table('/home/poirigui/51178_CMC_expmat.unfilt.data.txt.gz', compression='gzip', comment='#')
sc_df = sc_df.set_index(['Probe', 'Sequence', 'GeneSymbol', 'GeneName', 'GemmaId', 'NCBIid'])


# In[102]:


df = pd.read_table('/home/poirigui/7070_GSE37703_expmat.unfilt.data.txt.gz', compression='gzip', comment='#')
df = df.set_index(['Probe', 'Sequence', 'GeneSymbol', 'GeneName', 'GemmaId', 'NCBIid'])


# In[103]:


np.log10(sc_df.T.mean()).hist(bins=100)


# In[104]:


np.log10(df.T.mean()).hist(bins=100)


# In[119]:


import matplotlib.pyplot as plt


# In[125]:


np.pow(10, (-.5))


# In[126]:


np.log10(sc_df.T.var()).hist(bins=100)
plt.axvline(np.log10(0.5), c='r')


# In[ ]:





# In[122]:


np.log10(df.T.var()[df.T.var()>0.0001]).hist(bins=100)
plt.axvline(np.log10(0.05), c='r')


# In[ ]:




