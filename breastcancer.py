
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
from itertools import groupby
import datetime
from matplotlib import pyplot as plt




#show data
data=sc.read_10x_mtx('/home/xixi949592/breastcancer', cache=True)

df=data.transpose().to_df()

# Define column names and ranges
column_ranges = [
    ('BC01_', range(0, 3921)),
    ('BC02_', range(3921, 8797)),
    ('BC03_', range(8797, 13048)),
    ('BC04_', range(13048, 16606)),
    ('BC05_', range(16606, 21078)),
    ('BC06_', range(21078, 21279)),
    ('BC07_', range(21279, 22998)),
    ('BC08_', range(22998, 25425)),
    ('BC09_', range(25425, 30041)),
    ('BC10_', range(30041, 30943)),
    ('BC11_', range(30943, 35805)),
    ('BC12_', range(35805, 36167)),
    ('BC13_', range(36167, 39808)),
    ('BC14_', range(39808, 44024))
]

# Rename columns
for prefix, column_range in column_ranges:
    for i, j in enumerate(column_range):
        colname = f'{prefix}{i}'
        df.columns.values[j] = colname

# cid=data.obs

# data_new=datadf.mask(datadf<1,0).transpose()

# data_new

# filtering
sc.pp.filter_cells(data, min_genes=200)
sc.pp.filter_genes(data, min_cells=3)
data.var['mt'] = data.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.scatter(data, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(data, x='total_counts', y='n_genes_by_counts')

adata = data[data.obs.n_genes_by_counts < 6000, :]
adata = data[data.obs.pct_counts_mt < 25, :]

#adata
#data.obs['n_genes']

# normalization
sc.pp.normalize_total(adata, target_sum=None)

#log
sc.pp.log1p(adata)

# batch correction
sc.external.pp.mnn_correct(adata)


df_new=adata.to_df().transpose()

df_new

df_new.filter(regex='BC01',axis=1)

gene_list = ['BC01', 'BC02','BC03','BC04','BC05','BC06','BC07','BC08','BC09','BC10','BC11','BC12','BC13','BC14',]
BC = []
for i, gene_i in enumerate(gene_list):
    BC_i = np.repeat(gene_list[i], df_new.filter(regex=gene_i,axis=1).shape[1])
    BC.append(BC_i)


np.concatenate(BC)

g = np.concatenate(BC)


adata.obs['n_genes']=g

adata.obs['n_genes']

# PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata,color=['n_genes'],show=True)


# umap visualization
sc.pp.neighbors(adata, n_neighbors=60, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['n_genes'])

# TSNE visulization
sc.tl.tsne(adata)
sc.pl.tsne(adata,color=['n_genes'])

sc.pp.neighbors(adata, n_neighbors=60, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added='clusters', resolution=0.25)
sc.pl.umap(adata, color='clusters', add_outline=True, legend_loc='on data',
           legend_fontsize=12, legend_fontoutline=2,frameon=False,
           title='clustering of cells', palette='Set1')


# finding marker gene for cell types
sc.pp.neighbors(adata, n_neighbors=60,n_pcs=40)
sc.tl.leiden(adata,resolution=0.25)
sc.tl.rank_genes_groups(adata, 'leiden',method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False,fontsize=14)

# marker genes from paper
sc.tl.umap(adata)
sc.pl.umap(adata,color=['CD24','EPCAM','KRT19','PDGFRA','CD3D','LYZ','CLDN5','COL1A1','CD79A','LILRA4','CD1C','MS4A2'])

# marker genes we found
sc.tl.umap(adata)
sc.pl.umap(adata,color=['CD79A','KRT19','CD3E','COL1A2','MUCL1','TPSAB1','TYROBP','CLDN5','MZB1','LILRA4'])

marker_genes_dict = {'Cancer cell': ['KRT19','MUCL1'],
                     'T-cell': ['CD3E'],
                     'B-cell': ['MZB1', 'CD79A'],
                     'Fibroblasts': ['COL1A2'],
                     'Myeloid cell': ['TYROBP'],
                    'Mast cell':['TPSAB1'],
                    'Endothelial cell':['CLDN5'],
                    'Dendritic cell':['LILRA4']}

sc.pl.dotplot(adata, marker_genes_dict, 'clusters', dendrogram=True)

# create a dictionary to map cluster to annotation label
cluster2annotation = {
     '0': 'T-cell',
     '1': 'Fibroblast',
     '2': 'Myeloid cell',
     '3': 'Cancer cell',
     '4': 'Cancer cell',
     '5': 'Cancer cell',
     '6': 'T-cell',
     '7': 'Endothelial cell',
     '8': 'B-cell',
     '9': 'B-cell',
     '10':'Cancer cell',
     '11':'Mast cell',
     '12':'Cancer cell',
     '13':'Dendritic cell',
     '14':'Fibroblast',
     '15':'T-cell'}

# add a new `.obs` column called `cell type` by mapping clusters to annotation using pandas `map` function
adata.obs['cell type'] = adata.obs['clusters'].map(cluster2annotation).astype('category')

sc.pl.dotplot(adata, marker_genes_dict, 'cell type', dendrogram=True)

sc.pl.umap(adata, color='cell type',add_outline=True, legend_loc='on data',
           legend_fontsize=8, legend_fontoutline=2,frameon=False,
           title='clustering of cells', palette='Set1')

tumorcell=adata.obs[adata.obs['cell type']=='Cancer cell'].index
bc=adata.obs[adata.obs['cell type']=='B-cell'].index
tc=adata.obs[adata.obs['cell type']=='T-cell'].index
fb=adata.obs[adata.obs['cell type']=='Fibroblast'].index
edc=adata.obs[adata.obs['cell type']=='Endothelial cell'].index
ddc=adata.obs[adata.obs['cell type']=='Dendritic cell'].index
myc=adata.obs[adata.obs['cell type']=='Myeloid cell'].index
msc=adata.obs[adata.obs['cell type']=='Mast cell'].index

A_1=df_new.filter(regex=r'BC14',axis=1)
B_1=df_new.filter(regex=r'BC13',axis=1)
C_1=df_new.filter(regex=r'BC01|BC06|BC08|BC09',axis=1)
D_1=df_new.filter(regex=r'BC02|BC03|BC04|BC05|BC07|BC10|BC11|BC12',axis=1)

c4_col_series=pd.Series(np.array(msc))
D_1_col_series=D_1.columns.to_series()
s = pd.concat([c4_col_series,D_1_col_series])
res = s[s.duplicated(keep=False)].unique()
len(res)


fig,ax=plt.subplots()
cancertype=['Luminal A','Luminal B','HER2+','TNBC']
number=[[3675,10,75,176,157,0,110,13],[298,507,2111,50,50,18,392,215],[1231,359,3199,4249,831,25,1107,64],[7886,2087,8427,2505,1184,91,2768,54]]
value=np.array(number).astype(np.int)
year=[str(i) for i in cancertype]
v8=[i[0]+i[1]+i[2]+i[3]+i[4]+ i[5]+ i[6]+ i[7]for i in value]
v7=[i[1]+i[2]+i[3]+i[4]+ i[5]+ i[6]+ i[7] for i in value]
v6=[i[2]+i[3]+i[4]+ i[5]+ i[6]+ i[7] for i in value]
v5=[i[3]+i[4]+ i[5]+ i[6]+ i[7] for i in value]
v4=[i[4]+i[5]+ i[6]+ i[7]for i in value]
v3=[i[5]+ i[6]+ i[7]for i in value]
v2=[i[6]+ i[7]for i in value]
v1=[i[7]for i in value]

ax.barh(year,v8,color="thistle",height=0.5)
ax.barh(year,v7,color="plum",height=0.5)
ax.barh(year,v6,color="violet",height=0.5)
ax.barh(year,v5,color="orchid",height=0.5)
ax.barh(year,v4,color="mediumorchid",height=0.5)
ax.barh(year,v3,color="darkorchid",height=0.5)
ax.barh(year,v2,color="purple",height=0.5)
ax.barh(year,v1,color="indigo",height=0.5)

# marker_genes = ['CD24','CD74','KRT19','CD3E','EGFL7','COL1A1','SSR4','TYROBP','HLA-B','TPSAB1','IGFBP7','DCN']

# pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)

# new_cluster_names = ['Endothelial cells1 and B cells', 'Enteric neurons',
#     'Endothelial cells(inflammation)', 'Lymphocytes and Tumor1', 'Epithelial Cells and Ductal Cells',
#     'highly related to mestasis cells', 'Endothelial cells2', 'Olfactory epithelial cells and Megakaryocytes',
#     'Endothelial cells and B cells','Adipocytes', 'Luminal epithelial cells/Tumor2', 'Luminal epithelial cells/Tumor3','Fibroblasts']
# data_trans.rename_categories('leiden', new_cluster_names)

# sc.pl.tsne(adata, color='leiden', legend_loc='right margin', title='', frameon=True)

adata.obs['leiden']

# Initialize dictionaries to store indices for different cell types
tcell_indices = {}
bcell_indices = {}
cancer_indices = {}

# Define cluster numbers for each cell type
clusters = {
    'tcell': ['0', '6', '15'],
    'bcell': ['8', '9'],
    'cancer': ['3', '4', '5', '12', '10']
}

# Iterate over each cell type
for cell_type, cluster_numbers in clusters.items():
    # Initialize a list to store indices for the current cell type
    cell_indices = []
    # Iterate over cluster numbers for the current cell type
    for cluster_number in cluster_numbers:
        # Extract indices based on cluster number
        cluster_index = adata.obs[adata.obs['leiden'] == cluster_number].index
        # Append indices to the list
        cell_indices.append(cluster_index)
    # Concatenate indices and store in respective dictionary
    if cell_type == 'tcell':
        tcell_indices[cell_type] = np.concatenate(cell_indices)
    elif cell_type == 'bcell':
        bcell_indices[cell_type] = np.concatenate(cell_indices)
    elif cell_type == 'cancer':
        cancer_indices[cell_type] = np.concatenate(cell_indices)

# Extract data based on indices
tcell = df_new[np.concatenate(list(tcell_indices.values()))]
bcell = df_new[np.concatenate(list(bcell_indices.values()))]
cancer = df_new[np.concatenate(list(cancer_indices.values()))]

# Print or use the extracted data
print("T Cell Data:", tcell)
print("B Cell Data:", bcell)
print("Cancer Cell Data:", cancer)



B = cancer.filter(regex=r'BC13', axis=1)

# Extracting data for fibroblast cells
fibroblast_indices = np.concatenate([adata.obs[adata.obs['leiden'] == cluster].index 
                                      for cluster in ['1', '14']])
fibroblast = df_new[fibroblast_indices]

# Extracting data for mast cells
mast_indices = adata.obs[adata.obs['leiden'] == '11'].index
mast = df_new[mast_indices]

# Extracting data for dendritic cells
dc_indices = adata.obs[adata.obs['leiden'] == '13'].index
dc = df_new[dc_indices]

# Extracting data for endothelial cells
ec_indices = adata.obs[adata.obs['leiden'] == '7'].index
ec = df_new[ec_indices]

# Extracting data for myeloid cells
my_indices = adata.obs[adata.obs['leiden'] == '2'].index
my = df_new[my_indices]




# Exporting data to CSV files
cancer.to_csv('/home/xixi949592/only_cancer_cells_new_new.txt', sep='\t')
bcell.to_csv('/home/xixi949592/only_B_cells.txt', sep='\t')
tcell.to_csv('/home/xixi949592/only_T_cells.txt', sep='\t')
fibroblast.to_csv('/home/xixi949592/only_fibroblasts_new_new.txt', sep='\t')
Mast.to_csv('/home/xixi949592/only_mastcell.txt', sep='\t')
DC.to_csv('/home/xixi949592/only_dendritic.txt', sep='\t')
EC.to_csv('/home/xixi949592/only_endothelial_cell.txt', sep='\t')
MY.to_csv('/home/xixi949592/only_myeloid_cell.txt', sep='\t')



# scanpy read in only cancer cell data
cancerdata=sc.read_text('/home/xixi949592/only_cancer_cells_new_new.txt').transpose()

# subtypes in cancer cells
sc.tl.pca(cancerdata, svd_solver='arpack')
sc.pp.neighbors(cancerdata, n_neighbors=60, n_pcs=40)
sc.tl.umap(cancerdata)
sc.tl.leiden(cancerdata, key_added='clusters', resolution=0.15)

sc.pl.umap(cancerdata, color='clusters', add_outline=True, legend_loc='on data',
           legend_fontsize=12, legend_fontoutline=2,frameon=False,
           title='clustering of cells', palette='Set1')

cancerdata.obs['leiden']
cdf=cancerdata.to_df().transpose()


patientid = ['BC01', 'BC02','BC03','BC04','BC05','BC06','BC07','BC08','BC09','BC10','BC11','BC12','BC13','BC14',]
id = []
for i, gene_i in enumerate(patientid):
    id_i = np.repeat(patientid[i], cdf.filter(regex=gene_i,axis=1).shape[1])
    id.append(id_i)

np.concatenate(id)

n = np.concatenate(id)
cancerdata.obs['n_genes']=n



sc.tl.pca(cancerdata, svd_solver='arpack')
sc.pp.neighbors(cancerdata, n_neighbors=60, n_pcs=40)
sc.tl.umap(cancerdata)
sc.tl.leiden(cancerdata,resolution=0.15)

sc.pl.umap(cancerdata, color='n_genes', add_outline=True, legend_loc='on data',
           legend_fontsize=12, legend_fontoutline=2,frameon=False,
           title='clustering of cancer cells', palette='Set1')

sc.tl.tsne(cancerdata)
sc.pl.tsne(cancerdata,color=['n_genes'],legend_loc='on data')

# highly expressed genes in cancer cells
sc.pp.neighbors(cancerdata, n_neighbors=60,n_pcs=40)
sc.tl.leiden(cancerdata,resolution=0.15)
sc.tl.rank_genes_groups(cancerdata, 'leiden',method='t-test')
sc.pl.rank_genes_groups(cancerdata, n_genes=5, sharey=False,fontsize=14)

sc.tl.umap(cancerdata)
sc.pl.umap(cancerdata,color=['MUCL1','PKIB','IGF2','S100A11','XIST'])

marker_genes_cancer = {'C1_MUCL1': ['MUCL1'],
                     'C3_PKIB': ['PKIB'],
                     'C2_S100A11': ['S100A11'],
                      'C4_XIST':['XIST'], 'C5_IGF2':['IGF2']}

sc.pl.dotplot(cancerdata, marker_genes_cancer, 'clusters', dendrogram=True)

cancer_annotation = {
     '0': 'C3_PKIB',
     '1': 'C1_MUCL1',
     '2': 'C2_S100A11',
     '3': 'C2_S100A11',
     '4': 'C5_IGF2',
     '5': 'C4_XIST',
     '6': 'C2_S100A11',
     '7':'C2_S100A11'}

cancerdata.obs['cell type'] = cancerdata.obs['clusters'].map(cancer_annotation).astype('category')

sc.pl.umap(cancerdata, color='cell type', legend_loc='on data',
           frameon=False, legend_fontsize=10)

cancerdata.obs['cell type']

C1=cancerdata.obs[cancerdata.obs['cell type']=='C1_MUCL1'].index
C2=cancerdata.obs[cancerdata.obs['cell type']=='C2_S100A11'].index
C3=cancerdata.obs[cancerdata.obs['cell type']=='C3_PKIB'].index
C4=cancerdata.obs[cancerdata.obs['cell type']=='C4_XIST'].index
C5=cancerdata.obs[cancerdata.obs['cell type']=='C5_IGF2'].index


from matplotlib import pyplot as plt
fig,ax=plt.subplots()
cancertype=['Luminal A','Luminal B','HER2+','TNBC']
number=[[2,13,3500,160,0],[1,266,0,31,0],[6,1015,0,210,0],[2712,3264,0,636,1274]]
value=np.array(number).astype(np.int)
year=[str(i) for i in cancertype]
v1=[i[0]+i[1]+i[2]+i[3]+i[4] for i in value]
v2=[i[1]+i[2]+i[3]+i[4] for i in value]
v3=[i[2]+i[3]+i[4] for i in value]
v4=[i[3]+i[4] for i in value]
v5=[i[4] for i in value]
ax.barh(year,v1,color="violet",height=0.5)
ax.barh(year,v2,color="pink",height=0.5)
ax.barh(year,v3,color="skyblue",height=0.5)
ax.barh(year,v4,color="lightgreen",height=0.5)
ax.barh(year,v5,color="lavender",height=0.5)

value

# getting marker genes in cancer cells


# Define gene indices
gene_indices = [
    'MALAT1', 'COX6C', 'TMSB10', 'PKIB', 'ASAH1', 'LY6E', 'S100A6', 'TFF3', 'GATA3', 'MUCL1', 
    'CRABP1', 'CALML5', 'FDCSP', 'XIST', 'KRT19', 'ESR1', 'ERBB2', 'NEAT1', 'S100A11', 'XBP1', 
    'PTPRC', 'GSTP1', 'ITM2B', 'PFN1', 'IGF2'
]

# Create a list of gene indices
gene_index_lists = [cdf[cdf.index == gene].index for gene in gene_indices]

# Concatenate gene indices into a single array
all_gene_indices = np.concatenate(gene_index_lists)

# Extract data for selected genes
cancer_new = cdf.loc[all_gene_indices]


def find_duplicates(col_series, A_col_series):
    s = pd.concat([col_series, A_col_series])
    return len(s[s.duplicated(keep=False)].unique())

# Find duplicates for A
A = cancer_new.filter(regex=r'BC14', axis=1)
c1_col_series = pd.Series(np.array(C1))
A_col_series = A.columns.to_series()
c1_duplicates_A = find_duplicates(c1_col_series, A_col_series)

c2_col_series = pd.Series(np.array(C2))
c2_duplicates_A = find_duplicates(c2_col_series, A_col_series)

c3_col_series = pd.Series(np.array(C3))
c3_duplicates_A = find_duplicates(c3_col_series, A_col_series)

c4_col_series = pd.Series(np.array(C4))
c4_duplicates_A = find_duplicates(c4_col_series, A_col_series)

c5_col_series = pd.Series(np.array(C5))
c5_duplicates_A = find_duplicates(c5_col_series, A_col_series)

# Find duplicates for B
B = cancer_new.filter(regex=r'BC13', axis=1)
B_col_series = B.columns.to_series()
c1_duplicates_B = find_duplicates(c1_col_series, B_col_series)
c2_duplicates_B = find_duplicates(c2_col_series, B_col_series)
c3_duplicates_B = find_duplicates(c3_col_series, B_col_series)
c4_duplicates_B = find_duplicates(c4_col_series, B_col_series)
c5_duplicates_B = find_duplicates(c5_col_series, B_col_series)

# Find duplicates for C
C = cancer_new.filter(regex=r'BC01|BC06|BC08|BC09', axis=1)
C_col_series = C.columns.to_series()
c1_duplicates_C = find_duplicates(c1_col_series, C_col_series)
c2_duplicates_C = find_duplicates(c2_col_series, C_col_series)
c3_duplicates_C = find_duplicates(c3_col_series, C_col_series)
c4_duplicates_C = find_duplicates(c4_col_series, C_col_series)
c5_duplicates_C = find_duplicates(c5_col_series, C_col_series)

# Find duplicates for D
D = cancer_new.filter(regex=r'BC02|BC03|BC04|BC05|BC07|BC10|BC11|BC12', axis=1)
D_col_series = D.columns.to_series()
c1_duplicates_D = find_duplicates(c1_col_series, D_col_series)
c2_duplicates_D = find_duplicates(c2_col_series, D_col_series)
c3_duplicates_D = find_duplicates(c3_col_series, D_col_series)
c4_duplicates_D = find_duplicates(c4_col_series, D_col_series)
c5_duplicates_D = find_duplicates(c5_col_series, D_col_series)

# Concatenate the dataframes
heatmap = pd.concat([A, B, C, D], axis=1)

ax = plt.subplots(figsize=(15, 8))
ax=sns.heatmap(heatmap)

heatmap.transpose()

heatmap.index

genes = heatmap.transpose().index
types = ['Luminal A',]*int(A.shape[1]) + ['Luminal B',]*int(B.shape[1]) + ['HER2']*int(C.shape[1]) + ['TNBC',]*int(D.shape[1])
tuples = list(zip(types, genes))
index = pd.MultiIndex.from_tuples(tuples, names=['first', 'second'])
d = np.array(heatmap.transpose())
df = pd.DataFrame(d, index=index, columns =[np.array(heatmap.index)])

pd.DataFrame(types)



def add_line(ax, xpos, ypos):
    line = plt.Line2D([ypos, ypos+ .2], [xpos, xpos], color='black', transform=ax.transAxes)
    line.set_clip_on(False)
    ax.add_line(line)

def label_len(my_index,level):
    labels = my_index.get_level_values(level)
    return [(k, sum(1 for i in g)) for k,g in groupby(labels)]

def label_group_bar_table(ax, df):
    xpos = -.2
    scale = 1./df.index.size
    for level in range(df.index.nlevels):
        pos = df.index.size
        for label, rpos in label_len(df.index,level):
            add_line(ax, pos*scale, xpos)
            pos -= rpos
            lypos = (pos + .5 * rpos)*scale
            ax.text(xpos+.1, lypos, label, ha='center', transform=ax.transAxes)
        add_line(ax, pos*scale , xpos)
        xpos -= .2

# df = test_table()
df
fig = plt.figure(figsize = (10, 10))
ax = fig.add_subplot(111)
sns.heatmap(df)

#Below 3 lines remove default labels
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')

label_group_bar_table(ax, df)
fig.subplots_adjust(bottom=.1*df.index.nlevels)
plt.show()

df

row_colors= pd.DataFrame(['r']*int(A.shape[1]) + ['b']*int(B.shape[1]) + ['g']*int(C.shape[1]) + ['y']*int(D.shape[1]),index=heatmap.transpose().index)
# type(row_colors[0])

row_colors

sns.clustermap(heatmap.transpose(),
                   row_colors = row_colors,row_cluster= False)

# B cells
Bcell=sc.read_text('/home/xixi949592/only_B_cells.txt').transpose()

sc.tl.pca(Bcell, svd_solver='arpack')
sc.pp.neighbors(Bcell, n_neighbors=60, n_pcs=40)
sc.tl.umap(Bcell)
sc.tl.leiden(Bcell, key_added='clusters', resolution=0.15)

sc.pl.umap(Bcell, color='clusters', add_outline=True, legend_loc='on data',
           legend_fontsize=12, legend_fontoutline=2,frameon=False,
           title='clustering of cells', palette='Set1')

# marker genes in B cells
sc.pp.neighbors(Bcell, n_neighbors=60)
sc.tl.leiden(Bcell,resolution=1)
sc.tl.rank_genes_groups(Bcell, 'leiden',method='t-test')
sc.pl.rank_genes_groups(Bcell, n_genes=5, sharey=False,fontsize=14)

sc.tl.umap(Bcell)
sc.pl.umap(Bcell,color=['MZB1','CD37','HBB'])

marker_genes_bcell = {'C1_MZB1': ['MZB1'],
                     'C2_CD37': ['CD37'],
                     'C3_HBB': ['HBB']}

sc.pl.dotplot(Bcell, marker_genes_bcell, 'clusters', dendrogram=True)

bcell_annotation = {
     '0': 'C1_MZB1',
     '1': 'C2_CD37',
     '2': 'C2_CD37',
     '3':'C3_HBB'}

Bcell.obs['cell type'] = Bcell.obs['clusters'].map(bcell_annotation).astype('category')

sc.pl.umap(Bcell, color='cell type', legend_loc='on data',
           frameon=False, legend_fontsize=10)

bdf=Bcell.to_df().transpose()

# Define marker genes for B cells
marker_genes = ['FTH1', 'CD74', 'MZB1', 'CD37', 'IGKC', 'IGLC2', 'IGLL5', 'FKBP11', 'DERL3', 'SEC11C']

# Extract indices for marker genes
gene_indices = [bdf[bdf.index == gene].index for gene in marker_genes]

# Concatenate indices and extract corresponding data from bdf
bcell_new = bdf.loc[np.concatenate(gene_indices)]


A1=bcell_new.filter(regex=r'BC14',axis=1)
B1=bcell_new.filter(regex=r'BC13',axis=1)
C11=bcell_new.filter(regex=r'BC01|BC06|BC08|BC09',axis=1)
D1=bcell_new.filter(regex=r'BC02|BC03|BC04|BC05|BC07|BC10|BC11|BC12',axis=1)
heatmap1=pd.concat([A1,B1,C11,D1], axis=1)
row_colors_1= pd.DataFrame(['r']*int(A1.shape[1]) + ['b']*int(B1.shape[1]) + ['g']*int(C11.shape[1]) + ['y']*int(D1.shape[1]),index=heatmap1.transpose().index)
sns.clustermap(heatmap1.transpose(),
                   row_colors = row_colors_1,row_cluster= False)

C1=Bcell.obs[Bcell.obs['cell type']=='C1_MZB1'].index
C2=Bcell.obs[Bcell.obs['cell type']=='C2_CD37'].index
C3=Bcell.obs[Bcell.obs['cell type']=='C3_HBB'].index

c3_col_series=pd.Series(np.array(C3))
D1_col_series=D1.columns.to_series()
s = pd.concat([c3_col_series,D1_col_series])
res = s[s.duplicated(keep=False)].unique()
len(res)

from matplotlib import pyplot as plt
fig,ax=plt.subplots()
cancertype=['Luminal A','Luminal B','HER2+','TNBC']
number=[[2,8,0],[253,239,15],[194,164,1],[1003,1076,8]]
value=np.array(number).astype(np.int)
year=[str(i) for i in cancertype]
v1=[i[0]+i[1]+i[2] for i in value]
v2=[i[1]+i[2] for i in value]
v3=[i[2] for i in value]
ax.barh(year,v1,color="violet",height=0.5)
ax.barh(year,v2,color="pink",height=0.5)
ax.barh(year,v3,color="skyblue",height=0.5)

# T cells
Tcell=sc.read_text('/home/xixi949592/only_T_cells.txt').transpose()

sc.tl.pca(Tcell, svd_solver='arpack')
sc.pp.neighbors(Tcell, n_neighbors=60, n_pcs=40)
sc.tl.umap(Tcell)
sc.tl.leiden(Tcell, key_added='clusters', resolution=0.15)
sc.pl.umap(Tcell, color='clusters', add_outline=True, legend_loc='on data',
           legend_fontsize=12, legend_fontoutline=2,frameon=False,
           title='clustering of cells', palette='Set1')

# marker genes in T cells
sc.pp.neighbors(Tcell, n_neighbors=60,n_pcs=40)
sc.tl.leiden(Tcell,resolution=1)
sc.tl.rank_genes_groups(Tcell, 'leiden',method='t-test')
sc.pl.rank_genes_groups(Tcell, n_genes=5, sharey=False,fontsize=14)

sc.tl.umap(Tcell)
sc.pl.umap(Tcell,color=['LTB','CCL5','PTPRCAP','CXCL13'])

marker_genes_tcell = {'C1_LTB': ['LTB'],
                     'C3_PTPRCAP': ['PTPRCAP'],
                     'C2_CCL5': ['CCL5'],
                     'C4_CXCL13':['CXCL13']}

sc.pl.dotplot(Tcell, marker_genes_tcell, 'clusters', dendrogram=True)

tcell_annotation = {
     '1': 'C2_CCL5',
     '2': 'C3_PTPRCAP',
     '3': 'C4_CXCL13',
     '0':'C1_LTB'}

Tcell.obs['cell type'] = Tcell.obs['clusters'].map(tcell_annotation).astype('category')

sc.pl.umap(Tcell, color='cell type', legend_loc='on data',
           frameon=False, legend_fontsize=10)

tdf=Tcell.to_df().transpose()

# Define marker genes for T cells
marker_genes = ['LTB', 'CCL5', 'PTPRCAP', 'CXCL13', 'FOXP3', 'CD7', 'GZMK', 'NKG7', 'TNFRSF4', 'CST7', 'CD8A', 'CTSW', 'GNLY', 'GZMB']

# Extract indices for marker genes
gene_indices = [tdf[tdf.index == gene].index for gene in marker_genes]

# Concatenate indices and extract corresponding data from tdf
tcell_new = tdf.loc[np.concatenate(gene_indices)]

A2=tcell_new.filter(regex=r'BC14',axis=1)
B2=tcell_new.filter(regex=r'BC13',axis=1)
C22=tcell_new.filter(regex=r'BC01|BC06|BC08|BC09',axis=1)
D2=tcell_new.filter(regex=r'BC02|BC03|BC04|BC05|BC07|BC10|BC11|BC12',axis=1)
heatmap2=pd.concat([A2,B2,C22,D2], axis=1)
row_colors_2= pd.DataFrame(['r']*int(A2.shape[1]) + ['b']*int(B2.shape[1]) + ['g']*int(C22.shape[1]) + ['y']*int(D2.shape[1]),index=heatmap2.transpose().index)
sns.clustermap(heatmap2.transpose(),
                   row_colors = row_colors_2,row_cluster= False)

C1=Tcell.obs[Tcell.obs['cell type']=='C1_LTB'].index
C2=Tcell.obs[Tcell.obs['cell type']=='C2_CCL5'].index
C3=Tcell.obs[Tcell.obs['cell type']=='C3_PTPRCAP'].index
C4=Tcell.obs[Tcell.obs['cell type']=='C4_CXCL13'].index

c4_col_series=pd.Series(np.array(C4))
D2_col_series=D2.columns.to_series()
s = pd.concat([c4_col_series,D2_col_series])
res = s[s.duplicated(keep=False)].unique()
len(res)

from matplotlib import pyplot as plt
fig,ax=plt.subplots()
cancertype=['Luminal A','Luminal B','HER2+','TNBC']
number=[[44,31,0,0],[1206,878,23,4],[1475,1614,48,50],[3171,2501,2404,351]]
value=np.array(number).astype(np.int)
year=[str(i) for i in cancertype]
v1=[i[0]+i[1]+i[2]+i[3] for i in value]
v2=[i[1]+i[2]+i[3] for i in value]
v3=[i[2]+i[3] for i in value]
v4=[i[3] for i in value]
ax.barh(year,v1,color="violet",height=0.5)
ax.barh(year,v2,color="pink",height=0.5)
ax.barh(year,v3,color="skyblue",height=0.5)
ax.barh(year,v4,color="lightgreen",height=0.5)

# Fibroblasts
Fibro=sc.read_text('/home/xixi949592/only_fibroblasts_new_new.txt').transpose()

sc.tl.pca(Fibro, svd_solver='arpack')
sc.pp.neighbors(Fibro, n_neighbors=60, n_pcs=40)
sc.tl.umap(Fibro)
sc.tl.leiden(Fibro, key_added='clusters', resolution=0.1)
sc.pl.umap(Fibro, color='clusters', add_outline=True, legend_loc='on data',
           legend_fontsize=12, legend_fontoutline=2,frameon=False,
           title='clustering of cells', palette='Set1')

# marker genes in fibroblasts
sc.pp.neighbors(Fibro, n_neighbors=60,n_pcs=40)
sc.tl.leiden(Fibro,resolution=0.1)
sc.tl.rank_genes_groups(Fibro, 'leiden',method='t-test')
sc.pl.rank_genes_groups(Fibro, n_genes=5, sharey=False,fontsize=14)

sc.tl.umap(Fibro)
sc.pl.umap(Fibro,color=['DCN','IGFBP7','EEF1D','NEAT1'])

marker_genes_fibro = {'C1_DCN': ['DCN'],
                     'C3_EEF1D': ['EEF1D'],
                     'C2_IGFBP7': ['IGFBP7'],
                     'C4_NEAT1':['NEAT1']}

sc.pl.dotplot(Fibro, marker_genes_fibro, 'clusters', dendrogram=True)

fibro_annotation = {
     '1': 'C3_EEF1D',
     '2': 'C2_IGFBP7',
     '3': 'C4_NEAT1',
     '0':'C1_DCN'}
Fibro.obs['cell type'] = Fibro.obs['clusters'].map(fibro_annotation).astype('category')

sc.pl.umap(Fibro, color='cell type', legend_loc='on data',
           frameon=False, legend_fontsize=10)

fdf=Fibro.to_df().transpose()

# Define marker genes for fibroblast cells
marker_genes = ['DCN', 'EEF1D', 'NEAT1', 'PDGFRB', 'LUM', 'CTSK', 'C1S', 'AEBP1', 'COL1A2', 'IGFBP7', 'FN1', 'COL18A1']

# Extract indices for marker genes
gene_indices = [fdf[fdf.index == gene].index for gene in marker_genes]

# Concatenate indices and extract corresponding data from fdf
fibro_new = fdf.loc[np.concatenate(gene_indices)]   

A5=fibro_new.filter(regex=r'BC14',axis=1)
B5=fibro_new.filter(regex=r'BC13',axis=1)
C5=fibro_new.filter(regex=r'BC01|BC06|BC08|BC09',axis=1)
D5=fibro_new.filter(regex=r'BC02|BC03|BC04|BC05|BC07|BC10|BC11|BC12',axis=1)
heatmap5=pd.concat([A5,B5,C5,D5], axis=1)
row_colors_5= pd.DataFrame(['r']*int(A5.shape[1]) + ['b']*int(B5.shape[1]) + ['g']*int(C5.shape[1]) + ['y']*int(D5.shape[1]),index=heatmap5.transpose().index)
sns.clustermap(heatmap5.transpose(),
                   row_colors = row_colors_5,row_cluster= False)

C1=Fibro.obs[Fibro.obs['cell type']=='C1_DCN'].index
C2=Fibro.obs[Fibro.obs['cell type']=='C2_IGFBP7'].index
C3=Fibro.obs[Fibro.obs['cell type']=='C3_EEF1D'].index
C4=Fibro.obs[Fibro.obs['cell type']=='C4_NEAT1'].index

c4_col_series=pd.Series(np.array(C4))
D5_col_series=D5.columns.to_series()
s = pd.concat([c4_col_series,D5_col_series])
res = s[s.duplicated(keep=False)].unique()
len(res)

from matplotlib import pyplot as plt
fig,ax=plt.subplots()
cancertype=['Luminal A','Luminal B','HER2+','TNBC']
number=[[69,66,27,14],[24,16,9,1],[2043,449,1461,377],[1083,437,675,310]]
value=np.array(number).astype(np.int)
year=[str(i) for i in cancertype]
v1=[i[0]+i[1]+i[2]+i[3] for i in value]
v2=[i[1]+i[2]+i[3] for i in value]
v3=[i[2]+i[3] for i in value]
v4=[i[3] for i in value]
ax.barh(year,v1,color="violet",height=0.5)
ax.barh(year,v2,color="pink",height=0.5)
ax.barh(year,v3,color="skyblue",height=0.5)
ax.barh(year,v4,color="lightgreen",height=0.5)



# endothelial cells
endo=sc.read_text('/home/xixi949592/only_endothelial_cell.txt').transpose()

endo

sc.tl.pca(endo, svd_solver='arpack')
sc.pp.neighbors(endo, n_neighbors=60, n_pcs=40)
sc.tl.umap(endo)
sc.tl.leiden(endo, key_added='clusters', resolution=0.3)
sc.pl.umap(endo, color='clusters', add_outline=True, legend_loc='on data',
           legend_fontsize=12, legend_fontoutline=2,frameon=False,
           title='clustering of cells', palette='Set1')

# marker genes in dendtritic cells
sc.pp.neighbors(endo, n_neighbors=60,n_pcs=40)
sc.tl.leiden(endo,resolution=0.3)
sc.tl.rank_genes_groups(endo, 'leiden',method='t-test')
sc.pl.rank_genes_groups(endo, n_genes=5, sharey=False,fontsize=14)

sc.tl.umap(endo)
sc.pl.umap(endo,color=['ACKR1','ACTB','BTNL9','PODXL'])

marker_genes_den = {'C1_ACKR1': ['ACKR1'],
                     'C2_ACTB': ['ACTB'],
                     'C3_BTNL9': ['BTNL9'],
                     'C4_PODXL':['PODXL']}

sc.pl.dotplot(endo, marker_genes_den, 'clusters', dendrogram=True)

endo_annotation = {
     '0': 'C1_ACKR1',
     '1': 'C2_ACTB',
     '2': 'C3_BTNL9',
     '3':'C4_PODXL'}
endo.obs['cell type'] = endo.obs['clusters'].map(endo_annotation).astype('category')

sc.pl.umap(endo, color='cell type', legend_loc='on data',
           frameon=False, legend_fontsize=10)

ddf=endo.to_df().transpose()

# Define marker genes for endothelial cells
marker_genes = ['VWF', 'CD36', 'ACKR1', 'COL4A1', 'COL4A2', 'HSPG2', 'ACTB', 'BTNL9', 'PODXL']

# Extract indices for marker genes
gene_indices = [ddf[ddf.index == gene].index for gene in marker_genes]

# Concatenate indices and extract corresponding data from ddf
endo_new = ddf.loc[np.concatenate(gene_indices)]

A6=endo_new.filter(regex=r'BC14',axis=1)
B6=endo_new.filter(regex=r'BC13',axis=1)
C6=endo_new.filter(regex=r'BC01|BC06|BC08|BC09',axis=1)
D6=endo_new.filter(regex=r'BC02|BC03|BC04|BC05|BC07|BC10|BC11|BC12',axis=1)
heatmap6=pd.concat([A6,B6,C6,D6], axis=1)
row_colors_6= pd.DataFrame(['r']*int(A6.shape[1]) + ['b']*int(B6.shape[1]) + ['g']*int(C6.shape[1]) + ['y']*int(D6.shape[1]),index=heatmap6.transpose().index)
sns.clustermap(heatmap6.transpose(),
                   row_colors = row_colors_6,row_cluster= False)

C1=endo.obs[endo.obs['cell type']=='C1_ACKR1'].index
C2=endo.obs[endo.obs['cell type']=='C2_ACTB'].index
C3=endo.obs[endo.obs['cell type']=='C3_BTNL9'].index
C4=endo.obs[endo.obs['cell type']=='C4_PODXL'].index

c4_col_series=pd.Series(np.array(C4))
D6_col_series=D6.columns.to_series()
s = pd.concat([c4_col_series,D6_col_series])
res = s[s.duplicated(keep=False)].unique()
len(res)

from matplotlib import pyplot as plt
fig,ax=plt.subplots()
cancertype=['Luminal A','Luminal B','HER2+','TNBC']
number=[[14,9,23,111],[9,10,7,24],[208,220,282,121],[442,338,233,171]]
value=np.array(number).astype(np.int)
year=[str(i) for i in cancertype]
v1=[i[0]+i[1]+i[2]+i[3] for i in value]
v2=[i[1]+i[2]+i[3] for i in value]
v3=[i[2]+i[3] for i in value]
v4=[i[3] for i in value]
ax.barh(year,v1,color="violet",height=0.5)
ax.barh(year,v2,color="pink",height=0.5)
ax.barh(year,v3,color="skyblue",height=0.5)
ax.barh(year,v4,color="lightgreen",height=0.5)





# Myeloid cells
mye=sc.read_text('/home/xixi949592/only_myeloid_cell.txt').transpose()

mye

sc.tl.pca(mye, svd_solver='arpack')
sc.pp.neighbors(mye, n_neighbors=60, n_pcs=40)
sc.tl.umap(mye)
sc.tl.leiden(mye, key_added='clusters', resolution=0.3)
sc.pl.umap(mye, color='clusters', add_outline=True, legend_loc='on data',
           legend_fontsize=12, legend_fontoutline=2,frameon=False,
           title='clustering of cells', palette='Set1')

# marker genes in myeloid cells
sc.pp.neighbors(mye, n_neighbors=60,n_pcs=40)
sc.tl.leiden(mye,resolution=0.3)
sc.tl.rank_genes_groups(mye, 'leiden',method='t-test')
sc.pl.rank_genes_groups(mye, n_genes=5, sharey=False,fontsize=14)

sc.tl.umap(mye)
sc.pl.umap(mye,color=['CTSB','LAMP3','S100A8','VCAN','CD3E','LY6E'])

marker_genes_mye = {'C1_CD3E': ['CD3E'],
                     'C3_CTSB': ['CTSB'],
                     'C4_S100A8': ['S100A8'],
                     'C2_LY6E':['LY6E'],
                     'C5_LAMP3':'LAMP3'}

mye_annotation = {
     '0': 'C1_CD3E',
     '1': 'C2_LY6E',
     '2': 'C3_CTSB',
     '3':'C3_CTSB',
     '4':'C4_S100A8',
     '5':'C5_LAMP3'}
mye.obs['cell type'] = mye.obs['clusters'].map(mye_annotation).astype('category')

sc.pl.umap(mye, color='cell type', legend_loc='on data',
           frameon=False, legend_fontsize=10)

mdf=mye.to_df().transpose()

# Define marker genes for myeloid cells
marker_genes = ['CTSB', 'LY6E', 'S100A8', 'CD3E', 'LAMP3', 'CST3', 'VCAN', 'FOS', 'TYROBP']

# Extract indices for marker genes
gene_indices = [mdf[mdf.index == gene].index for gene in marker_genes]

# Concatenate indices and extract corresponding data from mdf
mye_new = mdf.loc[np.concatenate(gene_indices)]

A7=mye_new.filter(regex=r'BC14',axis=1)
B7=mye_new.filter(regex=r'BC13',axis=1)
C7=mye_new.filter(regex=r'BC01|BC06|BC08|BC09',axis=1)
D7=mye_new.filter(regex=r'BC02|BC03|BC04|BC05|BC07|BC10|BC11|BC12',axis=1)
heatmap7=pd.concat([A7,B7,C7,D7], axis=1)
row_colors_7= pd.DataFrame(['r']*int(A7.shape[1]) + ['b']*int(B7.shape[1]) + ['g']*int(C7.shape[1]) + ['y']*int(D7.shape[1]),index=heatmap7.transpose().index)
sns.clustermap(heatmap7.transpose(),
                   row_colors = row_colors_7,row_cluster= False)

C1=mye.obs[mye.obs['cell type']=='C1_CD3E'].index
C2=mye.obs[mye.obs['cell type']=='C2_LY6E'].index
C3=mye.obs[mye.obs['cell type']=='C3_CTSB'].index
C4=mye.obs[mye.obs['cell type']=='C4_S100A8'].index
C5=mye.obs[mye.obs['cell type']=='C5_LAMP3'].index

c5_col_series=pd.Series(np.array(C5))
D7_col_series=D7.columns.to_series()
s = pd.concat([c5_col_series,D7_col_series])
res = s[s.duplicated(keep=False)].unique()
len(res)

from matplotlib import pyplot as plt
fig,ax=plt.subplots()
cancertype=['Luminal A','Luminal B','HER2+','TNBC']
number=[[1,1,108,0,0],[316,1,40,18,17],[277,87,698,20,25],[613,1101,706,310,38]]
value=np.array(number).astype(np.int)
year=[str(i) for i in cancertype]
v1=[i[0]+i[1]+i[2]+i[3]+i[4] for i in value]
v2=[i[1]+i[2]+i[3]+i[4] for i in value]
v3=[i[2]+i[3]+i[4] for i in value]
v4=[i[3]+i[4] for i in value]
v5=[i[4] for i in value]
ax.barh(year,v1,color="violet",height=0.5)
ax.barh(year,v2,color="pink",height=0.5)
ax.barh(year,v3,color="skyblue",height=0.5)
ax.barh(year,v4,color="lightgreen",height=0.5)
ax.barh(year,v5,color="lavender",height=0.5)



