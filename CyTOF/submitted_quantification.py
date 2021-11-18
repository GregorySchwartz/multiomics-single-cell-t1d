#!/usr/bin/env python
# coding: utf-8

# In[54]:


import pandas as pd
import numpy as np
import scipy.stats as stats


# In[18]:


from scipy.stats import chisquare


# ## 1. Quantify prop of disease state in each cell type

# In[20]:


import matplotlib.pyplot as plt
import pandas as pd
from dependencies.helper_functions import *
data_dir = '/mnt/data2/aanchal/data/cyTOF_T1D/'
dataname = 'cyTOF_T1D'


# In[21]:


classifier='ELM'


suffix= 'TUNE(no-1)_cofactor=200_thlist=1_unknownth=75'
sig_suffix ='_usingPanel'
annospat_cytof0 = pd.read_csv(data_dir+'labels_ssc/trte_labels_ssc_'+classifier+suffix+'.csv')

suffix='_AnnoSpatFINAL'
annospat_cytof = pd.read_csv(data_dir+'labels_ssc/trte_labels_ssc_'+classifier+suffix+'.csv')


# In[22]:


raw_data_file,signatures,roi_name_col,disease_status_col,p1,p2,toDrop,suffix,data_dir,subfolder_name,gt = get_data(dataname, sig_suffix)


# In[23]:


cell_status = raw_data_file['Status'].reset_index().rename(columns = {'index':'item'})
cell_status


# In[24]:


cell_label_status = annospat_cytof0 .merge(cell_status,on ='item' )
cell_label_status_immunesubtypes =annospat_cytof.merge(cell_status,on ='item' )


# In[25]:


cell_label_status.label.unique()


# In[26]:


cell_label_status_immunesubtypes.label.unique()


# In[27]:


cell_label_status.Status.unique()


# #### Donor grps percentage per cell type

# In[28]:


count = cell_label_status.groupby(['label','Status']).count()#.apply(lambda x: 100 * x / float(x.sum()))#['Status'].count()#.unstack().plot.bar(stacked=True, figsize=(17, 10))
count_immunesubtypes = cell_label_status_immunesubtypes.groupby(['label','Status']).count()#.apply(lambda x: 100 * x / float(x.sum()))#['Status'].count()#.unstack().plot.bar(stacked=True, figsize=(17, 10))


# In[29]:


count.groupby(level=0).apply(lambda x:100 * x / float(x.sum())).unstack().plot.bar(stacked=True, figsize=(15, 3), width=0.65,color = ['#FFA500','#008000','#4B0082'])#plt.show()
plt.legend(loc='best', bbox_to_anchor=(1, 0.5))
plt.savefig("figures/stackedBar_DonorGroupPercentagePerCT_"+dataname+".pdf", bbox_inches='tight')


# In[30]:


"figures/stackedBar_DonorGroupPercentagePerCT_"+dataname+".pdf"


# In[31]:


count_immunesubtypes.groupby(level=0).apply(lambda x:100 * x / float(x.sum())).unstack().plot.bar(stacked=True, figsize=(15, 3), width=0.65,color = ['#FFA500','#008000','#4B0082'])#plt.show()
plt.legend(loc='best', bbox_to_anchor=(1, 0.5))
plt.savefig("figures/stackedBar_DonorGroupPercentagePerCT_"+dataname+"_immunesubtypes.pdf", bbox_inches='tight')


# In[32]:


"figures/stackedBar_DonorGroupPercentagePerCT_"+dataname+"_immunesubtypes.pdf"


# ## 2. Counting fraction of disease state in beta, immune, ductal and hybrid

# #### cells in a branch subsetted using: grep -e ,98, -e /81/ -e /87/ tmc_annospatoverlay.log > beta_cells.log
# 

# In[33]:


celltype='Hybrid'########Dont use "branch subsets" for annotated cell types#########'Ductal' #'Immune'#'Beta'#'Hybrid'#' Immune' 'Ductal' 'Hybrid'
ct_cells_ids = pd.read_csv('/mnt/data2/aanchal/analysis/NatMetaReview_aries/cytof/'+celltype+'_cells.log', delimiter = "\t", header=None)
ct_cells_ids = ct_cells_ids[0].str.extract(pat='(^.*?,)')[0].str.replace(",","")
#prolif_cell_ids = prolif_cells_ids[0].str.extract(pat='(^.*?-)')[0].str.replace("-","")
subset = cell_label_status[cell_label_status['item'].isin(ct_cells_ids)]


#OR

celltype='Ductal_Cells'#'Beta_Cells'#'Immune cells'
subset = cell_label_status[cell_label_status['label']==celltype]
#arr=subset.Status.value_counts().reindex(['T1D','GAD+','GAD-']).rename(index={'GAD+':'AAB+','GAD-':'Control'})
subset.Status.value_counts()


# In[34]:



subset.Status.value_counts()


# In[35]:


arr=subset.Status.value_counts()


# In[36]:


totalcells=np.sum(arr)
statistic,p =chisquare(arr) #https://www.statology.org/chi-square-goodness-of-fit-test-python/, https://www.statology.org/chi-square-p-value-calculator/ statistic=chisquare score

if (p==0):
    p='<1e-6'
else:
    p='='+str(p)


# In[37]:


pd.DataFrame(arr, index=['T1D','AAB+','Control']).transpose()


# In[38]:


arr


# In[39]:


tab = pd.DataFrame(arr, index=['T1D','AAB+','Control']).transpose().loc['Status']#pd.crosstab(df['CTQ-tool'],df['opinion']).loc['d']


# In[40]:


tab


# In[41]:


ax = tab.plot.pie(startangle=90,colors = ['#4B0082','#FFA500','#008000'])
ax.set_ylabel('')
ax.set_title('cyTOF\n'+celltype+'\n(n='+str(totalcells)+', p'+p+')', fontweight='bold')
plt.savefig("figures/pie_DonorGroupDivisionIn_"+celltype+"_"+dataname+".pdf", bbox_inches='tight')


# In[42]:


"figures/pie_DonorGroupDivisionIn_"+celltype+"_"+dataname+".pdf"


# ## Quantify the percentage of hybrid cells in each disease condition (T1D/AAB+/Control)

# In[43]:


import seaborn as sns; sns.set_theme(); sns.set(color_codes=True)
from statannot import add_stat_annotation
sns.set(style="whitegrid")


# In[44]:


cell_label_status_temp = cell_label_status.copy()
cell_label_status_temp_immunesubtypes =  cell_label_status_immunesubtypes.copy()


# In[45]:


celltype='Hybrid'#'Immune' 'Ductal' 'Hybrid'
ct_cells_ids = pd.read_csv('/mnt/data2/aanchal/analysis/NatMetaReview_aries/cytof/'+celltype+'_cells.log', delimiter = "\t", header=None)
ct_cells_ids = ct_cells_ids[0].str.extract(pat='(^.*?,)')[0].str.replace(",","")


# In[46]:


#cell_label_status[
cell_label_status_temp.loc[cell_label_status_temp['item'].isin(ct_cells_ids ) ,'label']='Hybrid_Cells'
cell_label_status_temp_immunesubtypes.loc[cell_label_status_temp['item'].isin(ct_cells_ids ) ,'label']='Hybrid_Cells'
#]


# In[47]:


cell_label_status_temp['Donor'] =cell_label_status_temp.item.str.partition('-')[2]#.extract(pat='(-^.*?)')
cell_label_status_temp_immunesubtypes['Donor'] =cell_label_status_temp_immunesubtypes.item.str.partition('-')[2]#.extract(pat='(-^.*?)')


# In[48]:


cell_label_status_temp


# In[49]:


df_unstacked = cell_label_status_temp.groupby(['Donor','Status','label']).size().unstack().reset_index()
df_unstacked_immunesubtypes= cell_label_status_temp_immunesubtypes.groupby(['Donor','Status','label']).size().unstack().reset_index()


# In[50]:


df2plot = df_unstacked.groupby(['Donor','Status']).sum().apply(lambda x:100 * x / float(x.sum()), axis=1).reset_index()
#df2plot['log(1+Hybrid_Cells)']=np.log( 100+df2plot.Hybrid_Cells.astype('int') )
df2plot_immunesubtypes = df_unstacked_immunesubtypes.groupby(['Donor','Status']).sum().apply(lambda x:100 * x / float(x.sum()), axis=1).reset_index()


# In[ ]:





# In[51]:


ax=sns.boxplot(data=df2plot, x="Status", y="Hybrid_Cells", hue="Status")
#ax._legend.remove()
ax = sns.swarmplot(data=df2plot, x="Status", y="Hybrid_Cells", hue="Status", dodge=True, color=".25")
# add_stat_annotation(ax, data=df2plot, x="Status", y="Hybrid_Cells", hue="Status",
#                     box_pairs=[ (("T1D", "AAB+")) ],
#                    test='Kruskal', loc='inside' ,verbose=2)

plt.legend(loc='best', bbox_to_anchor=(1.3, 1))
plt.savefig("figures/boxplot_HybridPercentageOfTotal_"+dataname+".pdf")


# In[56]:


df2plot_for1wayANOVA = df2plot.pivot(index='Donor',columns='Status',values='Hybrid_Cells')
# https://www.reneshbedre.com/blog/anova.html

fvalue, pvalue = stats.f_oneway(df2plot_for1wayANOVA['T1D'].dropna().values, df2plot_for1wayANOVA['AAB+'].dropna().values, df2plot_for1wayANOVA['Control'].dropna().values)
print(fvalue, pvalue)
df2plot_for1wayANOVA


# In[57]:


df2plot.loc[df2plot['Status']=='Control'][['Donor','Hybrid_Cells']]
#cytOF


# In[58]:


df2plot


# In[59]:


ax=sns.boxplot(data=df2plot, x="Status", y="Immune cells", hue="Status")
#ax._legend.remove()
ax = sns.swarmplot(data=df2plot, x="Status", y="Immune cells", hue="Status", dodge=True, color=".25")
# add_stat_annotation(ax, data=df2plot, x="Status", y="Hybrid_Cells", hue="Status",
#                     box_pairs=[ (("T1D", "AAB+")) ],
#                    test='Kruskal', loc='inside' ,verbose=2)

plt.legend(loc='best', bbox_to_anchor=(1.3, 1))
plt.savefig("figures/boxplot_ImmunePercentageOfTotal_"+dataname+".pdf")


# In[60]:


df2plot.loc[df2plot['Status']=='Control'][['Donor','Immune cells']]


# In[ ]:





# In[ ]:





# ### S17-D avg of (Hybrid cell %age of total cells per donor) across all donors in disease state (T1D,AAB+,COntrol)

# In[61]:


import scipy #https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.stats.mstats.linregress.html


# In[62]:


df2plot2=df2plot[["Donor","Immune cells","Hybrid_Cells","Status"]]
df2plot2


# In[63]:


df2plot2_immunesubtypes=df2plot_immunesubtypes[["Donor",'Hybrid_Cells',"Status",'Helper_T_Cells','Killer_T_Cells', 'Memory_T_Cells', 'Myeloid_Cells','NK_Cells', 'Regulatory_T_Cells']]        
df2plot2_immunesubtypes


# In[64]:


#df2plot2.groupby(['Status']).mean().reindex(["Control","AAB+","T1D"])
a1 =df2plot2.groupby(['Status']).mean().reindex(["Control","AAB+","T1D"])
a1_immunesubtypes =df2plot2_immunesubtypes.groupby(['Status']).mean().reindex(["Control","AAB+","T1D"])
a2 =df2plot2.groupby(['Status']).median().reindex(["Control","AAB+","T1D"])
a2_immunesubtypes =df2plot2_immunesubtypes.groupby(['Status']).median().reindex(["Control","AAB+","T1D"])

ct2='Hybrid_Cells'
ct1='Immune cells'


# In[65]:


a=a1.copy()

x=a[ct1]
y=a[ct2]
import scipy #https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.stats.mstats.linregress.html
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
r_squared = r_value**2
plt.plot(x, y, 'o', label='original data')
plt.plot(x, intercept + slope*x, 'r', label='fitted line')
for idx, row in a.iterrows(): 
    plt.text(row['Immune cells'], row['Hybrid_Cells'], idx)
plt.text(np.max(a['Immune cells']), np.min(a['Hybrid_Cells']), "R^2 = "+str(np.round(r_squared,5) ) )
plt.savefig('figures/corr_meanHybridVsImmune%.pdf')


# In[66]:


a=a2.copy()

x=a[ct1]
y=a[ct2]
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
r_squared = r_value**2
plt.plot(x, y, 'o', label='original data')
plt.plot(x, intercept + slope*x, 'r', label='fitted line')
for idx, row in a.iterrows(): 
    plt.text(row['Immune cells'], row['Hybrid_Cells'], idx)
    
plt.text(np.max(a['Immune cells']), np.min(a['Hybrid_Cells']), "R^2 = "+str(np.round(r_squared,5) ) )
plt.savefig('figures/corr_medianHybridVsImmune%.pdf')


# In[67]:


for ct1 in ['Killer_T_Cells', 'Helper_T_Cells','Memory_T_Cells', 'Myeloid_Cells','NK_Cells', 'Regulatory_T_Cells']:
    print(ct1)
    
    a=a1_immunesubtypes.copy()
    x=a[ct1]
    y=a[ct2]
    import scipy #https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.stats.mstats.linregress.html
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
    r_squared = r_value**2
    plt.plot(x, y, 'o', label='original data')
    plt.plot(x, intercept + slope*x, 'r', label='fitted line')
    for idx, row in a.iterrows(): 
        plt.text(row[ct1], row[ct2], idx)
    plt.text(np.max(a[ct1]), np.min(a[ct1]), "R^2 = "+str(np.round(r_squared,5) ) )
    plt.title("Hybrid vs "+ct1)
    plt.savefig('figures/corr_meanHybridVs'+ct1+'%.pdf')
    plt.show()


# In[68]:


for ct1 in ['Killer_T_Cells', 'Helper_T_Cells','Memory_T_Cells', 'Myeloid_Cells','NK_Cells', 'Regulatory_T_Cells']:
    print(ct1)
    
    a=a2_immunesubtypes.copy()
    x=a[ct1]
    y=a[ct2]
    import scipy #https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.stats.mstats.linregress.html
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
    r_squared = r_value**2
    plt.plot(x, y, 'o', label='original data')
    plt.plot(x, intercept + slope*x, 'r', label='fitted line')
    for idx, row in a.iterrows(): 
        plt.text(row[ct1], row[ct2], idx)
    plt.text(np.max(a[ct1]), np.min(a[ct1]), "R^2 = "+str(np.round(r_squared,5) ) )
    plt.title("Hybrid vs "+ct1)
    plt.savefig('figures/corr_medianHybridVs'+ct1+'%.pdf')
    plt.show()


# In[ ]:





# In[ ]:




