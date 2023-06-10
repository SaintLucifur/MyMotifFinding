#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import matplotlib.pyplot as plt

# Load p-values from the first dataset
with open('sox2HomerResult.txt', 'r') as file:
    homer_pvalues = [float(line.strip()) for line in file.readlines()[1:]]

# Load p-values from the second dataset
with open('sox2MyMotifFindingResults.txt', 'r') as fn:
    mymotif_pvalues = [float(line.strip()) for line in fn.readlines()[1:]]

top_10_home_pvalues = sorted(homer_pvalues)[:10]
top_10_mymotif_pvalues = sorted(mymotif_pvalues)[:10]
    
x = range(1, len(top_10_home_pvalues) + 1)

plt.figure(figsize=(8, 6))
plt.semilogy(x, top_10_home_pvalues, 'b', label='homer result')
plt.semilogy(x, top_10_mymotif_pvalues, 'r', label='MyMotifFinding result')
plt.xlabel('Index')
plt.ylabel('P-values')
plt.title('Top 10 P-values Comparison for SOX2')
plt.legend()

# Adjust the y-axis range
plt.ylim([1e-300, 1e-10])

# Show the plot
plt.show()

