#!/usr/bin/env python
# coding: utf-8

# In[119]:


#Uncomment lineA to merge immune cell subtypes


# In[120]:


import pandas as pd
import numpy as np


# In[121]:


from scipy.stats import chisquare


# In[122]:


import matplotlib.pyplot as plt
dataname = 'IMC_T1D'
data_dir = '/mnt/data2/aanchal/data/'+dataname+'/'

suffix ='_FINAL'
suffix ='_withImmuneCelltypes_withNegMarkers_d-adaptive2_HPAP026Control_9'


# In[123]:


classifier='ELM'
ssc_labels = pd.read_csv(data_dir+'labels_ssc/trte_labels_ssc_'+classifier+suffix+'.csv')

# lineA below:
ssc_labels0 = ssc_labels.copy().replace({'Memory_T_Cells':'Immune_Cells' , 'Regulatory_T_Cells':'Immune_Cells', 'Killer_T_Cells':'Immune_Cells','Helper_T_Cells':'Immune_Cells', 'NK_Cells':'Immune_Cells','B_Cells':'Immune_Cells', 'Myeloid_Cells':'Immune_Cells'})

raw_data_file = pd.read_csv('/mnt/data2/aanchal/data/IMC_T1D/raw_data/mgDF.csv', index_col=0)
raw_data_file.loc[raw_data_file.TIFFfilename.str.startswith('HPAP-026'),'Status']='GAD-'

cell_status = raw_data_file['Status'].reset_index().rename(columns = {'index':'item'})

cell_label_status_immunesubtypes = ssc_labels.merge(cell_status,on ='item' )
cell_label_status = ssc_labels0.merge(cell_status,on ='item' )


# In[124]:


from dependencies.helper_functions import map_roi_names

cell_label_status_temp = cell_label_status.copy()


raw_data_df=pd.read_csv('/mnt/data2/aanchal/data/IMC_T1D/raw_data/mgDF.csv', index_col=0)
raw_data_df['Donor_Part_ROI']    = map_roi_names(raw_data_df[['TIFFfilename']])['Donor_Part_ROI']
raw_data_df['Donor']  = raw_data_df['Donor_Part_ROI'].str.extract(pat='([H][P][A][P]...)')
raw_data_df['Part']  = raw_data_df['Donor_Part_ROI'].str[8:12]

#cell_donor = raw_data_df[['Donor']].reset_index().rename(columns ={'index':'item'})
cell_donor_part = raw_data_df[['Donor','Part']].reset_index().rename(columns ={'index':'item'})
cell_label_status_part  = cell_label_status_temp.merge(cell_donor_part, left_on ='item', right_on ='item')
cell_label_status_part_immunesubtypes = cell_label_status_immunesubtypes.merge(cell_donor_part, left_on ='item', right_on ='item')


# In[125]:


cell_label_status_part.to_csv('/mnt/data2/aanchal/data/IMC_T1D/raw_data/item_label_Status_Donor_Part_HPAP026shifted2control.csv')
cell_label_status_part 


# In[126]:


cell_label_status_part.label.unique()


# ## Pie chart

# In[127]:



# #acinar
# celltype='Acinar'
# arr=[1115,294,20043]

# #beta
# celltype='Beta'
# arr=[471,13726,16396]

#proliferating
# celltype='Proliferating'
# arr=[4333,3614,1543]

#hybrid
celltype='HLA-DR+ductal'
#arr=[34983, 5712, 2634]
celltype='Hybrid_Cells'#'Immune' 'Ductal' 'Hybrid'
ct_cells_ids = pd.read_csv('/mnt/data2/aanchal/analysis/NatMetaReview_aries/imc/logs/'+celltype+'.log', delimiter = "\t", header=None)
ct_cells_ids = ct_cells_ids[0].str.extract(pat='(^.*?-)')[0].str.replace("-","")
cell_label_status.loc[cell_label_status['item'].isin(ct_cells_ids ) ,'label']=celltype


#celltype='Beta_Cells'#'Immune_Cells' #'Immune_Cells' # 'Ductal_Cells'
subset = cell_label_status[cell_label_status['label']==celltype]
arr=subset.Status.value_counts().reindex(['T1D','GAD+','GAD-']).rename(index={'GAD+':'AAB+','GAD-':'Control'})
arr


# In[128]:


cell_label_status.label.unique()#[cell_label_status['label']==celltype]


# In[129]:


totalcells=np.sum(arr)

statistic,p =chisquare(arr) #https://www.statology.org/chi-square-goodness-of-fit-test-python/, https://www.statology.org/chi-square-p-value-calculator/ statistic=chisquare score


# In[130]:


if (p==0):
    p='<1e-6'
else:
    p='='+str(p)


# In[131]:


pd.DataFrame(arr, index=['T1D','AAB+','Control']).transpose().loc['Status']


# In[ ]:





# In[132]:


tab = pd.DataFrame(arr, index=['T1D','AAB+','Control']).transpose().loc['Status']#pd.crosstab(df['CTQ-tool'],df['opinion']).loc['d']

ax = tab.plot.pie(startangle=90,colors = ['#4B0082','#FFA500','#008000'])
ax.set_ylabel('')
ax.set_title('IMC\n'+celltype+'\n(n='+str(totalcells)+', p'+p+')', fontweight='bold')
plt.savefig("figures/pie_DonorGroupDivisionIn_"+celltype+"_"+dataname+"_HPAP026shifted2control.pdf")


# In[133]:


"figures/pie_DonorGroupDivisionIn_"+celltype+"_"+dataname+"_HPAP026shifted2control.pdf"


# ### Which cell type is most proliferating ?

# In[134]:


#subset prolif cells from branch using >grep ,112, tmc_annospatoverlay.log > proliferating_cells.log
prolif_cells_ids = pd.read_csv('/mnt/data2/aanchal/analysis/NatMetaReview_aries/imc/logs/proliferating_cells.log', delimiter = "\t", header=None)

prolif_cell_ids = prolif_cells_ids[0].str.extract(pat='(^.*?-)')[0].str.replace("-","")


# In[135]:


dataname = 'IMC_T1D'
data_dir = '/mnt/data2/aanchal/data/'+dataname+'/'

suffix = '_withImmuneCelltypes_withNegMarkers_d-adaptive2'
suffix ='_FINAL'
suffix ='_withImmuneCelltypes_withNegMarkers_d-adaptive2_HPAP026Control_9'

classifier='ELM'
#ssc_labels0 = pd.read_csv(data_dir+'labels_ssc/trte_labels_ssc_'+classifier+suffix+'.csv')


# In[136]:


ssc_labels_prolif_cells = ssc_labels0[ssc_labels0['item'].isin(prolif_cell_ids )]


# In[137]:


ssc_labels0[ssc_labels0['item'].isin(prolif_cell_ids )]['label'].value_counts().plot.bar()


# In[138]:


ssc_labels[ssc_labels['item'].isin(prolif_cell_ids )]['label'].value_counts().plot.bar()


# ## Donor grps percentage per cell type

# In[139]:


count = cell_label_status.groupby(['label','Status']).count()#.apply(lambda x: 100 * x / float(x.sum()))#['Status'].count()#.unstack().plot.bar(stacked=True, figsize=(17, 10))
count_immunesubtypes = cell_label_status_immunesubtypes.groupby(['label','Status']).count()#.apply(lambda x: 100 * x / float(x.sum()))#['Status'].count()#.unstack().plot.bar(stacked=True, figsize=(17, 10))


# In[140]:


count.groupby(level=0).apply(lambda x:100 * x / float(x.sum())).unstack().plot.bar(stacked=True, figsize=(15, 3), width=0.65,color = ['#FFA500','#008000','#4B0082'])#plt.show()
plt.legend(loc='best', bbox_to_anchor=(1, 0.5))
plt.savefig("figures/stackedBar_DonorGroupPercentagePerCT_"+dataname+"_immuneSubtypesMreged_HPAP026shifted2control.pdf", bbox_inches='tight')


# In[141]:


"figures/stackedBar_DonorGroupPercentagePerCT_"+dataname+"_immuneSubtypesMreged_HPAP026shifted2control.pdf"


# In[142]:


count_immunesubtypes.groupby(level=0).apply(lambda x:100 * x / float(x.sum())).unstack().plot.bar(stacked=True, figsize=(15, 3), width=0.65,color = ['#FFA500','#008000','#4B0082'])#plt.show()
plt.legend(loc='best', bbox_to_anchor=(1, 0.5))
plt.savefig("figures/stackedBar_DonorGroupPercentagePerCT_"+dataname+"_HPAP026shifted2control.pdf", bbox_inches='tight')


# In[143]:


"figures/stackedBar_DonorGroupPercentagePerCT_"+dataname+"_HPAP026shifted2control.pdf"


# ## Quantify the percentage of hybrid cells in each disease condition (T1D/AAB+/Control)

# In[144]:


cell_label_status_part


# In[146]:


ct_cells_ids = pd.read_csv('/mnt/data2/aanchal/analysis/NatMetaReview_aries/imc/logs/Hybrid_Cells.log', delimiter = "\t", header=None)
ct_cells_ids = ct_cells_ids[0].str.extract(pat='(^.*?-)')[0].str.replace("-","")

cell_label_status_part.loc[cell_label_status_part['item'].isin(ct_cells_ids ) ,'label']='Hybrid_Cells'
cell_label_status_part_immunesubtypes.loc[cell_label_status_part['item'].isin(ct_cells_ids ) ,'label']='Hybrid_Cells'


# In[147]:


ct_cells_ids


# In[148]:


df_unstacked = cell_label_status_part.groupby(['Donor','Status','label']).size().unstack().reset_index()
df_unstacked_immunesubtypes= cell_label_status_part_immunesubtypes.groupby(['Donor','Status','label']).size().unstack().reset_index()


# In[149]:


df2plot = df_unstacked.groupby(['Donor','Status']).sum().apply(lambda x:100 * x / float(x.sum()), axis=1).reset_index()
df2plot_immunesubtypes = df_unstacked_immunesubtypes.groupby(['Donor','Status']).sum().apply(lambda x:100 * x / float(x.sum()), axis=1).reset_index()


# In[150]:


import seaborn as sns; sns.set_theme(); sns.set(color_codes=True)
from statannot import add_stat_annotation
sns.set(style="whitegrid")


# In[151]:


ax=sns.boxplot(data=df2plot, x="Status", y="Hybrid_Cells", hue="Status")
#ax._legend.remove()
ax = sns.swarmplot(data=df2plot, x="Status", y="Hybrid_Cells", hue="Status", dodge=True, color=".25")
# add_stat_annotation(ax, data=df2plot, x="Status", y="Hybrid_Cells", hue="Status",
#                     box_pairs=[ (("GAD-", "GAD+")) ],
#                    test='t-test_ind', loc='inside' ,verbose=2)

plt.legend(loc='best', bbox_to_anchor=(2.5, 1))
plt.savefig("figures/boxplot_HybridPercentageOfTotal_"+dataname+"_HPAP026shifted2control.pdf", bbox_inches='tight')


# In[152]:


"figures/boxplot_HybridPercentageOfTotal_"+dataname+"_HPAP026shifted2control.pdf"


# In[153]:


import scipy.stats as stats
df2plot_for1wayANOVA = df2plot.pivot(index='Donor',columns='Status',values='Hybrid_Cells')
#https://www.reneshbedre.com/blog/anova.html

fvalue, pvalue = stats.f_oneway(df2plot_for1wayANOVA['T1D'].dropna(), df2plot_for1wayANOVA['GAD+'].dropna(), df2plot_for1wayANOVA['GAD-'].dropna())
print(fvalue, pvalue)


# In[154]:


df2plot.loc[df2plot['Status']=='GAD-']
#IMC


# In[155]:


ax=sns.boxplot(data=df2plot, x="Status", y="Immune_Cells", hue="Status")
#ax._legend.remove()
ax = sns.swarmplot(data=df2plot, x="Status", y="Immune_Cells", hue="Status", dodge=True, color=".25")
# add_stat_annotation(ax, data=df2plot, x="Status", y="Hybrid_Cells", hue="Status",
#                     box_pairs=[ (("GAD-", "GAD+")) ],
#                    test='t-test_ind', loc='inside' ,verbose=2)

plt.legend(loc='best', bbox_to_anchor=(2.5, 1))
plt.savefig("figures/boxplot_ImmunePercentageOfTotal_"+dataname+"_HPAP026shifted2control.pdf", bbox_inches='tight')


# In[156]:


"figures/boxplot_ImmunePercentageOfTotal_"+dataname+"_HPAP026shifted2control.pdf"


# In[157]:


import scipy.stats as stats
df2plot_for1wayANOVA = df2plot.pivot(index='Donor',columns='Status',values='Immune_Cells')
#https://www.reneshbedre.com/blog/anova.html

fvalue, pvalue = stats.f_oneway(df2plot_for1wayANOVA['T1D'].dropna(), df2plot_for1wayANOVA['GAD+'].dropna(), df2plot_for1wayANOVA['GAD-'].dropna())
print(fvalue, pvalue)


# In[158]:


df2plot.loc[df2plot['Status']=='GAD+']


# ### S17-E- avg of (Hybrid cell %age of total cells per donor) across all donors in disease state (T1D,AAB+,COntrol)

# In[159]:


import scipy #https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.stats.mstats.linregress.html


# In[160]:


df2plot.columns


# In[161]:


df2plot_immunesubtypes


# In[162]:


df2plot


# In[163]:


df2plot2=df2plot[["Donor","Immune_Cells","Hybrid_Cells","Status"]]

#df2plot2_immunesubtypes=df2plot_immunesubtypes[["Donor",'Hybrid_Cells',"Status",'Helper_T_Cells','Killer_T_Cells', 'Memory_T_Cells', 'Myeloid_Cells','NK_Cells', 'Regulatory_T_Cells']]        
df2plot2


# In[164]:


a1 =df2plot2.groupby(['Status']).mean().reindex(["GAD-","GAD+","T1D"])
a1_immunesubtypes =df2plot2_immunesubtypes.groupby(['Status']).mean().reindex(["GAD-","GAD+","T1D"])

a2 =df2plot2.groupby(['Status']).median().reindex(["GAD-","GAD+","T1D"])
a2_immunesubtypes =df2plot2_immunesubtypes.groupby(['Status']).median().reindex(["GAD-","GAD+","T1D"])

ct2='Hybrid_Cells'
ct1='Immune_Cells'


# In[165]:


#for ct1 in ['Immune_Cells']:
print(ct1)

a=a1.copy()
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
plt.savefig('figures/corr_meanHybridVs'+ct1+'%_HPAP026shifted2control.pdf')
plt.show()


# In[166]:


#for ct1 in ['Immune_Cells']:
print(ct1)

a=a2.copy()
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
plt.savefig('figures/corr_medianHybridVs'+ct1+'%_HPAP026shifted2control.pdf')
plt.show()


# In[167]:


'figures/corr_HybridVs'+ct1+'_HPAP026shifted2control.pdf'


# In[168]:


a2_immunesubtypes


# In[169]:


a1_immunesubtypes


# In[170]:


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
    plt.savefig('figures/corr_medianHybridVs'+ct1+'%_HPAP026shifted2control.pdf')
    plt.show()


# In[171]:


'figures/corr_HybridVs'+ct1+'_HPAP026shifted2control.pdf'


# In[172]:


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
    plt.savefig('figures/corr_meanHybridVs'+ct1+'%_HPAP026shifted2control.pdf')
    plt.show()


# In[173]:


# a=a2_immunesubtypes.copy()

# x=a[ct1]
# y=a[ct2]
# slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
# r_squared = r_value**2
# plt.plot(x, y, 'o', label='original data')
# plt.plot(x, intercept + slope*x, 'r', label='fitted line')
# for idx, row in a.iterrows(): 
#     plt.text(row[ct1], row[ct2], idx)
    
# plt.text(np.max(a[ct1])-25, np.min(a[ct2]), "R^2 = "+str(np.round(r_squared,5) ) )


# In[ ]:





# ### S17-C- Plot Hrbrid cell %age of total cells per pancreatic region (fig S17C)

# In[56]:


cell_label_status_part


# In[57]:


df_unstacked2 = cell_label_status_part.groupby(['Part','label']).size().unstack().reset_index()
df_unstacked2


# In[58]:


df2plot


# In[59]:


df2plot = df_unstacked2.groupby(['Part']).sum().apply(lambda x:100 * x / float(x.sum()), axis=1)#.reset_index()
df2plot=df2plot[['Hybrid_Cells']]
df2plot


# In[62]:


df2plot.reindex(["Head", "Body", "Tail","None"]).plot.bar(color = ['#808080'], rot=45)
plt.legend(loc='best', bbox_to_anchor=(1, 0.5))
plt.ylabel("Percentage of Total Cells %")
plt.xlabel("")
plt.savefig("figures/Bar_HybridPercentageOfTotalPerPancreaticregion_"+dataname+"_HPAP026shifted2control.pdf", bbox_inches='tight')


# In[64]:


"figures/Bar_HybridPercentageOfTotalPerPancreaticregion_"+dataname+"_HPAP026shifted2control.pdf"


# In[ ]:





# In[65]:


#s17D,E--incomplete


# In[66]:


cell_label_status_part


# In[67]:


df_unstacked2 = cell_label_status_part.groupby(['Status','label']).size().unstack().reset_index()
df_unstacked2


# In[68]:


df2plot = df_unstacked2.groupby(['Status']).sum().apply(lambda x:100 * x / float(x.sum()), axis=1)#.reset_index()
df2plot=df2plot[['Hybrid_Cells']]


# In[69]:


df2plot


# In[ ]:





# ### Which CTs are around hybrid ?

# In[71]:


#read csv file labelling the neighbors
neigh_info = pd.read_csv('/mnt/data2/aanchal/analysis/NatMetaReview_aries/imc/output_gw/output_hpap_imc_5_cut/output_node_257_proximity_20/labels.csv')
neigh_info = neigh_info.set_index(keys='item')

cell_id_by_gw = pd.read_csv('/mnt/data2/aanchal/data/IMC_T1D/raw_data/imc_fromGW.csv')['item']
neigh_info = neigh_info.reindex(cell_id_by_gw).reset_index()


# In[72]:


#cell_label_status_temp.merge(neigh_info, left_on='item', right_on='item')
cell_label_status_part['Neighbor_label'] = neigh_info['label']
cell_label_status_part


# In[73]:


cell_label_status_part.Neighbor_label.unique()


# In[74]:


cell_label_status_part.columns


# In[75]:


#subset the "Neighbors"
base = cell_label_status_part.loc[cell_label_status_part['Neighbor_label'] =='Base']
neighbors = cell_label_status_part.loc[cell_label_status_part['Neighbor_label'] =='Neighbor']
neighbors_immunesubtypes = cell_label_status_part_immunesubtypes.loc[cell_label_status_part['Neighbor_label'] =='Neighbor']
neighbors
#ssc_labels0[ssc_labels0['item'].isin(prolif_cell_ids )]['label'].value_counts().plot.bar()


# In[76]:


#check the composition of neighbors
celltype = 'Neighbor'#'Base'
arr=neighbors.Status.value_counts().reindex(['T1D','GAD+','GAD-']).rename(index={'GAD+':'AAB+','GAD-':'Control'})

totalcells=np.sum(arr)
statistic,p =chisquare(arr)
if (p==0):
    p='<1e-6'
else:
    p='='+str(p)
tab = pd.DataFrame(arr, index=['T1D','AAB+','Control']).transpose().loc['Status']#pd.crosstab(df['CTQ-tool'],df['opinion']).loc['d']
ax = tab.plot.pie(startangle=90,colors = ['#4B0082','#FFA500','#008000'])
ax.set_ylabel('')
ax.set_title('IMC\n'+celltype+' cells\n(n='+str(totalcells)+', p'+p+')', fontweight='bold')


# In[97]:


neighbors['label'].value_counts().plot.bar()
plt.savefig("figures/bar_CTsAroundHybrid_"+dataname+"_immuneSubtypesMreged_HPAP026shifted2control.pdf", bbox_inches='tight')


# In[100]:


"figures/bar_CTsAroundHybrid_"+dataname+"_immuneSubtypesMreged_HPAP026shifted2control.pdf"


# In[99]:


neighbors_immunesubtypes['label'].value_counts().plot.bar()
plt.savefig("figures/bar_CTsAroundHybrid_"+dataname+"_HPAP026shifted2control.pdf", bbox_inches='tight')


# In[79]:


neighbors


# In[80]:


count = neighbors.groupby(['Neighbor_label','Status']).count()#.apply(lambda x: 100 * x / float(x.sum()))#['Status'].count()#.unstack().plot.bar(stacked=True, figsize=(17, 10))
count#.groupby(level='Neighbor_label').apply(lambda x:100 * x / float(x.sum())).unstack().plot.bar(stacked=True, figsize=(15, 3), width=0.65,color = ['#FFA500','#008000','#4B0082'])#plt.show()

