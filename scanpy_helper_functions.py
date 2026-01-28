### Helper Functions for Scanpy Workflows:
import pandas as pd
import scanpy as sc
import numpy as np
import scipy.sparse
import anndata as an
import harmonypy # type: ignore
import matplotlib as plt
import matplotlib.colors as mcolors
import random
from pybiomart import Server, Dataset # type: ignore
import sklearn

def check_adata_format(adatas, batches, data_sources, celltype_colnames, subtype_colnames, condition_colnames, age_colnames, sex_colnames, save_intermediate_path=None):
    '''
    Input:
        - adatas: list of adata objects
        - batches: list of batch column names in same order as adata objects
        - data_sources: list of data source names in same order as adata objects
        - celltype_colnames: list of cell type column names in same order as adata objects
        - subtype_colnames: list of cell type column names in same order as adata objects (set as NA if there is no subtype column available)
        - condition_colnames: list of condition/phenotype column names in the same order as adata objects
        - age_colnames: list of age column names in the same order as adata objects (set with values of NA if not available)
        - sex_colnames: list of sex column names in the same order as adata objects (set with values of NA if not available)
        - save_intermediate: each passed adata file will be saved to the following path (added later, helpful if joining adatas is running out of memory)
    '''
    return_adatas=[]
    for i, adata in enumerate(adatas):
        adata=adatas[i]
        batch=batches[i]
        data_source=data_sources[i]
        celltype_colname=celltype_colnames[i]
        subtype_colname=subtype_colnames[i]
        condition_colname=condition_colnames[i]
        age_colname=age_colnames[i]
        sex_colname=sex_colnames[i]

        print(data_source, 'trying to reset counts matrix to raw data...')
        if 'counts' in adata.layers:
            print(data_source, ' adata.layers["counts"] exists... resetting counts \n')
            adata.X=adata.layers['counts']
        elif hasattr(adata.raw, 'X'):
            print(data_source, ' adata.raw.X exists... resetting counts \n')
            adata.X=adata.raw.X
        elif adata.X.max()>20:
            print('counts matrix is already likely unprocessed. Max counts currently: ', adata.X.max(), "\n")
        else:
            print(data_source, ' no raw data available, continuing with counts as is \n')
        
        adata.var_names_make_unique()

        ######################################################### Filtering:
        print('standard preprocess...\n')
        print(data_source, 'pre-filtering steps shape: ', adata.shape, "\n")
        adata.var["mt"]=adata.var_names.str.startswith("MT-") #human mitochondrial genes specified
        adata.var["ribo"]=adata.var_names.str.startswith(("RPS", "RPL"))
        adata.var["hb"]=adata.var_names.str.contains("^HB[^(P)]")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
        
        sc.pp.filter_cells(adata, min_genes=100)
        sc.pp.filter_genes(adata, min_cells=3)
        print('number of cells with greater than 10% mitochondrial genes',len(adata[(adata.obs.pct_counts_mt > 10)]), "\n")
        adata=adata[(adata.obs.pct_counts_mt < 10)] #filtering out cells with greater than 10% mitochondrial genes
        
        #Ensuring that batch is category type
        if adata.obs[batch].dtype.name != 'category':
            adata.obs[batch]=adata.obs[batch].apply(lambda x: f'Batch{x}' if adata.obs[batch].dtype.name != 'category' else x).astype('category')
        ######################################################### Scrublet to Identify Doublets:
        print('starting scrublet analysis')
        print('pre-scrublet adata shape: ', adata.shape, "\n")
        if batch != None:
            try:
                sc.pp.scrublet(adata, batch_key=batch)
                predicted_doublet_idx=list(adata[adata.obs.predicted_doublet==True].obs.index)
                print(len(predicted_doublet_idx), " number of doublets predicted across the dataset ", data_source, "\n")
                if len(predicted_doublet_idx) > 0:
                    adata=adata[~adata.obs.index.isin(predicted_doublet_idx)] #filtering cells in predicted doublet index list out of adata object
                    print('doublets removed from', data_source, 'new shape: ', adata.shape, "\n")
                else:
                    print('no doublets to remove from', data_source, "\n")
            except Exception as e:
                print('scrublet with batch but without batch separation failed for ', data_source, "\n")
                print(f'Error: {str(e)}', "\n")
                batch_ids=adata.obs[batch].unique()
                predicted_doublet_idx=[]
                for batch in batch_ids:
                    print("scrubbing batch: ", batch)
                    adata_batch=adata[adata.obs["batch"]==batch]
                    sc.pp.scrublet(adata_batch) #performing scrublet on individual batches #may need to look at adjusting number of principal components later

                    predicted_doublet_idx.append(list(adata_batch[adata_batch.obs.predicted_doublet==True].obs.index)) #pulling a list of cell where predicted_doublet is true and adding this to the list of predicted doublet cells for the entire batch
    
                flat_predicted_doublet_idx=[cell_id for batch in predicted_doublet_idx for cell_id in batch]
                print(len(flat_predicted_doublet_idx), ": number of doublets predicted across the dataset ", data_source, "\n")     
                if len(flat_predicted_doublet_idx) > 0:
                    print('Done with scrublet, starting filter out of doublets \n')
                    adata=adata[~adata.obs.index.isin(flat_predicted_doublet_idx)] #filtering cells in predicted doublet index list out of adata object
                    print('done filtering of doublets from adata ', data_source, "\n")
                else: #no cells to filter out so continue
                    print(data_source, ' has no predicted doublets. Adata shape: ', adata.shape, "\n")
        else: #batch being None
            sc.pp.scrublet(adata)
            predicted_doublet_idx=list(adata[adata.obs.predicted_doublet==True].obs.index)
            print(len(predicted_doublet_idx), " number of doublets predicted across the dataset ", data_source, '\n')
            adata=adata[~adata.obs.index.isin(predicted_doublet_idx)] #filtering cells in predicted doublet index list out of adata object
            print('done filtering of doublets from adata ', data_source)
        print('adata shape after scrublet: ', adata.shape, '\n')
        ################################################# Normalizing Data:
        if adata.X.max()<20:
            print(data_source,' already logged and normalized counts \n')
        else:
            print(adata, 'normalizing and log transforming... \n')
            adata.layers["counts"]=adata.X.copy()
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
        
        ############################################# Feature Selection:
        try: 
            adata.var['highly_variable']
            print(i, ' already has highly variable genes \n')
        except:
            if batch!='None':
                    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=batch)
            else:
                sc.pp.highly_variable_genes(adata, n_top_genes=2000)
                
        ############################################# Standardizing Obs Columns:
        if 'datasource' not in adata.obs.columns:
            adata.obs['datasource']=np.full(adata.obs.shape[0], data_source)

        if celltype_colname != 'Celltype':
            adata.obs.rename(columns={celltype_colname: 'Celltype'}, inplace=True)

        if subtype_colname != None:
            if subtype_colname != 'Subtype':
                adata.obs.rename(columns={subtype_colname: 'Subtype'}, inplace=True)
        else:
            adata.obs['Subtype']=np.full(adata.obs.shape[0], 'NA')

        if condition_colname != 'Condition':
            adata.obs.rename(columns={condition_colname: 'Condition'}, inplace=True)
        #Checking number of conditions (MS has subtypes of MS in the Condition col)
        try:
            if adata.obs['Condition'].nunique() > 3:
                adata.obs.rename(columns={'Condition': 'Original_Condition'}, inplace=True)
        
                adata.obs['Condition']=adata.obs['Original_Condition'].map(
                    lambda x: 'Control' if x=='CTR' else ('MS' if x in ['PPMS', 'RRMS', 'SPMS'] else x)
                    )
            else:
                print('\n Condition column okay with conditions:', list(adata.obs.Condition.unique()), "\n")
        except Exception as e:
            print(f"\n Error updating Condition column: {e} \n")

        if age_colname != None:
            if age_colname != 'Age':
                adata.obs.rename(columns={age_colname: 'Age'}, inplace=True)
                adata.obs.Age=adata.obs.Age.round().astype('Int64')
        else:
            adata.obs['Age']=pd.Series([pd.NA] * adata.obs.shape[0], dtype='Int64')
        
        standardized_sex_encoding={'male': 'M', 'Male': 'M', 'MALE': 'M', 'm': 'M', 'M': 'M', 1: 'M',
                                     'female': 'F', 'Female': 'F', 'FEMALE': 'F', 'f': 'F', 'F': 'F', 0: 'F', 'XX': 'F', 'XY': 'M', 'XYY': 'M'}
        common_categories=['M', 'F', 'NA'] 
        if sex_colname != None:
            if sex_colname != 'Sex':
                mapped_sex=adata.obs[sex_colname].map(standardized_sex_encoding)
                mapped_sex=mapped_sex.astype("object").fillna('NA')
                adata.obs['Sex']=pd.Categorical(mapped_sex, categories=common_categories)
        else:
            adata.obs['Sex']=pd.Categorical(['NA'] * adata.shape[0], categories=common_categories)

        #asserting categorical for Sex
        adata.obs['Sex']=adata.obs['Sex'].astype(pd.CategoricalDtype(categories=common_categories))

        ############################################# PCA + Integration:
        if batch!=None:
            try:
                sc.tl.pca(adata)
                sc.external.pp.scanorama_integrate(adata, key=batch)
            except Exception as e:
                error_msg=str(e)
                print(f"{data_source}: Scanorama error occurred: {error_msg} \n")
                if "non-contiguous batches" in error_msg.lower():
                    print(f"{data_source}: Attempting to fix non-contiguous batches by sorting...")
                    try:
                        adata=adata[np.argsort(adata.obs[batch].values)].copy()
                        sc.external.pp.scanorama_integrate(adata, key=batch)
                        print(f"{data_source}: Scanorama successfully rerun with contiguous batches.\n")
                    except Exception as inner_e:
                        print(f"{data_source}: Retry failed after sorting: {inner_e} \n")
        else:
            print('No scanorama since no batch given \n')
        #writing to intermediate file if save_intermediate_path:
        if save_intermediate_path!=None:
            cur_save_intermediate_path=save_intermediate_path[i]
            adata.write(cur_save_intermediate_path)
        return_adatas.append(adata)

    return return_adatas

def assign_celltype_class(celltype):
    '''
    Intended to be used to assign a dataset's original celltype definition into the AD dataset's broader
    celltype definitions.

    Apply with code such as this: 
    adata.obs['Celltype_Class']=adata.obs.apply(
    lambda row: sh.assign_celltype_class(row['Celltype']), axis=1
)
    '''
    if isinstance(celltype, str) and (celltype.lower().startswith("exc") or celltype.lower().startswith('ex-') or celltype.lower().startswith('exn')):
        return "Excitatory neurons"
    elif isinstance(celltype, str) and (celltype.lower().startswith("inh") or celltype.lower().startswith('in-') or celltype.lower().startswith('inn')):
        return "Inhibitory neurons"
    elif isinstance(celltype, str) and celltype.lower().startswith("ast"):
        return "Astrocytes"
    elif isinstance(celltype, str) and (celltype.lower().startswith("mic") or celltype.lower().startswith("mg") or celltype.lower()=='t cells'):
        return "Microglia"
    elif isinstance(celltype, str) and (celltype.lower().startswith("oli") or celltype.lower().startswith("odc")):
        return "Oligodendrocytes"
    elif isinstance(celltype, str) and (celltype.lower().startswith('vas') or celltype.lower().startswith('endo') or celltype.lower().startswith('peri')):
        return "Vascular Cells"
    elif isinstance(celltype, str) and celltype.lower().startswith("opc"):
        return "OPCs + COPs"
    else:
        return celltype
    
def assign_celltype_fullnames(celltype):
    '''
    This function should be used to take the Celltype from AD and create a new column with full names (more friendly to plotting)
    Apply with code such as this: 
    adata.obs['Celltype_Class']=adata.obs.apply(
    lambda row: sh.assign_celltype_fullnames(row['Celltype']), axis=1
    )
    '''
    if isinstance(celltype, str) and celltype.lower().startswith("exc"):
        return "Excitatory neurons"
    elif isinstance(celltype, str) and celltype.lower().startswith("inh"):
        return "Inhibitory neurons"
    elif isinstance(celltype, str) and celltype.lower().startswith("ast"):
        return "Astrocytes"
    elif isinstance(celltype, str) and celltype.lower().startswith("mic"):
        return "Microglia"
    elif isinstance(celltype, str) and celltype.lower().startswith("oli"):
        return "Oligodendrocytes"
    elif isinstance(celltype, str) and celltype.lower().startswith('vas'):
        return "Vascular Cells"
    elif isinstance(celltype, str) and celltype.lower().startswith("opc"):
        return "OPCs + COPs"
    else:
        return celltype
    
def create_umaps(adata, adata_name, colnames, date, sep_by=None, sep_by_colnames=None, extra_info=None, point_size=2):
    '''
    Input:
        - adata: assumption is that adata has been run through check_adata_format() and join_adatas() first
        - adata_name: dataset name/identifier (string)
        - colnames: list of column names to be colored by
        - date: date as a string
        - sep_by: .obs column to separate adatas into and plot separately. Chosen column must be categorical.
        - sep_by_colnames: columns to be colored by in UMAP
        - extra_info: additional information to be added into the saved png name
    '''
    if sep_by!=None:
        print('trying to split', sep_by,' for sep by umap')
        split_categories=list(adata.obs[sep_by].unique())
        for i in split_categories:
            cur_adata=adata[adata.obs[sep_by]==i]
            try: 
                for col in sep_by_colnames:
                    if extra_info!=None:
                        sep_save_str=adata_name+"_"+i+"_"+col+"_"+date+"_"+extra_info+'.png'
                    else:
                        sep_save_str=adata_name+"_"+i+"_"+col+"_"+date+'.png'
                    i_title=i.replace('_', ' ')+': '+col
                    sc.pl.umap(cur_adata, color=col, size=point_size, title=i_title, save=sep_save_str)
            except Exception as e:
                print('could not save umap(s) for', i, "due to...")
                print(f"{e}")

    for colname in colnames:
        print('starting umap for...', colname)
        if extra_info!=None:
            save_str=adata_name+"_"+colname+"_"+date+"_"+extra_info+'.png'
        else:
            save_str=adata_name+"_"+colname+"_"+date+'.png'
        sc.pl.umap(adata, color=colname, size=point_size, save=save_str)
        
    return f"Done saving umaps for {adata_name}"

def mtx_to_adata(counts_path, gene_names_path, metadata_path, genename_col=None, barcodes_col=None, data_name=None, directory_path=None, raw_counts_path=None, spatial_path=None):
    '''
    Input:
        - counts: full path to matrix of preprocessed counts (or if only raw counts are available then raw counts)
        - gene_names: full path to csv of gene names to be set as varnames (must be same length and order as matrix)
        - metadata: full path to csv of obs information (length should correspond to number of cells)
        - barcodes_col: column name for the barcode/row names in the metadata (if None then no row names will be set)
        - data_name: string with name of data for file saving
        - directory_path: path to location for the new file to be stored
        - raw_counts: path to matrix of raw counts (use if both preprocessed and raw counts are available)
    '''
    adata=sc.read_mtx(counts_path)
    if metadata_path.endswith('.csv') or metadata_path.endswith('.txt'):
        metadata=pd.read_csv(metadata_path)
    elif metadata_path.endswith('.tsv'):
        metadata=pd.read_csv(metadata_path, sep='\t')
    if gene_names_path.endswith('.csv') or gene_names_path.endswith('.txt'):
        genes=pd.read_csv(gene_names_path)
    elif gene_names_path.endswith('.tsv'):
        genes=pd.read_csv(gene_names_path, sep='\t')

    if len(genes)==adata.shape[1] and len(metadata)==adata.shape[0]:
        print('processed counts matrix orientation is correct') #no need to transpose
    elif len(genes)==adata.shape[0] and len(metadata)==adata.shape[1]:
        print('processed counts matrix is flipped, will need to transpose')
        adata=adata.T
    
    if spatial_path != None:
        spatial=pd.read_csv(spatial_path)
        if spatial.shape[0] != adata.shape[0]:
            print('Warning: The number of rows in spatial data does not match the number of cells in the count matrix.')
        else:
            print('Spatial data successfully loaded')
            # Save spatial coordinates to adata.obsm
            adata.obsm['spatial']=spatial.values
    #checking metadata columns for mixed types
    for col in metadata.columns:
        if metadata[col].apply(type).nunique() > 1:
            print(f"Column '{col}' has mixed types. Converting to string.")
            metadata[col]=metadata[col].astype(str)
    adata.obs=metadata
    if genename_col!=None:
        adata.var_names=genes[genename_col]
    else:
        adata.var_names=genes.iloc[:, 0].values
    if barcodes_col!=None:
        adata.obs_names=metadata[barcodes_col]

    if raw_counts_path !=None:
        raw=sc.read_mtx(raw_counts_path)
        if len(genes)==raw.shape[1] and len(metadata)==raw.shape[0]:
            print('raw counts matrix orientation is correct') #no need to transpose
        elif len(genes)==raw.shape[0] and len(metadata)==raw.shape[1]:
            print('raw counts matrix is flipped, will need to transpose')
            raw=raw.T
        else:
            print('shape of raw counts matrix and genes/metadata does not match')
            exit()
        adata.layers['raw counts']=raw

    if data_name != None and directory_path!=None:
        adata.write(directory_path+data_name+'.h5ad')
    
    print('Created Adata...', adata)
    print('Max counts value: ', adata.X.max())
    return(adata)

def add_gene_names_to_adata(adata):
    #Dictionary mapping Ensembl IDs to gene names from pybiomart import Server, Dataset
    server=Server(host='http://www.ensembl.org')
    dataset=Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    ensembl_df=dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
    ensembl_df=ensembl_df[ensembl_df['Gene stable ID'].isin(list(adata.var_names))==True].dropna()
    df=pd.Series(ensembl_df['Gene name'].values, index=ensembl_df['Gene stable ID'].values).to_dict()

    #Add 'gene_name' to adata
    gene_names=[]
    valid_ensembl_ids=[]
    for ensembl_id in adata.var_names:
        if ensembl_id in df:
            gene_names.append(df[ensembl_id])
            valid_ensembl_ids.append(ensembl_id)
            
    adata.var['gene_name']=pd.Series(gene_names, index=valid_ensembl_ids)
    origin_genes_num=adata.shape[1]
    #Remove the genes that are not in the Ensembl df from the adata
    adata=adata[:, adata.var_names.isin(valid_ensembl_ids)]
    print("Number of genes removed due to no mapping gene id: ", origin_genes_num-adata.shape[1])
    
    return adata

def plot_celltype_bar(adata, x_col, save_string, colorby=None, percent_true=False, title=None, x_label=None, y_label=None):
    if x_label is None:
        x_label=x_col.replace("_", " ").title()
    if y_label is None:
        y_label="Percentage of Cells" if percent_true else "Number of Cells"
    if title is None:
        title=f"{'Proportion' if percent_true else 'Number'} of Cells per {x_label}"

    if colorby:
        counts=adata.obs.groupby(x_col, observed=False)[colorby].value_counts()
        counts_unstacked=counts.unstack(fill_value=0)
        if percent_true:
            data=counts_unstacked.div(counts_unstacked.sum(axis=1), axis=0) * 100
        else:
            data=counts_unstacked
        ax=data.plot(kind='bar', stacked=True, figsize=(10, 6))
        legend_title=colorby.replace("_", " ").title()
        ax.legend(title=legend_title, bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        counts=adata.obs[x_col].value_counts().sort_index()
        if percent_true:
            data=counts / counts.sum() * 100
        else:
            data=counts
        ax=data.plot(kind='bar', color='skyblue', figsize=(8, 5))

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    plt.tight_layout()
    plt.figsave(save_string)

def set_celltype_colors(adata=None, celltype_class_col=None, celltype_subtype_col=None, broad=False, condition_col=None, datasource_col=None):
    '''
    Used to make colors consistent for subtypes and celltypes across datasets
    '''
    broad_colors={
        'Excitatory neurons':'#006d2c', #deep green
        'OPCs + COPs': '#8c510a', #brown
        'Vascular Cells': '#bebada', #purple
        'Astrocytes': '#FF7518', #orange
        'Inhibitory neurons': '#377eb8', #blue
        'Oligodendrocytes': '#f781bf', #soft pink
        'Microglia': '#e41a1c' #red
        }
    subtype_colors={
        #Excitatory neurons (green)
        'Exc L2-3 CBLN2 LINC02306': "#00441b",'Exc L3-4 RORB CUX2': "#006d2c",'Exc L3-5 RORB PLCH1': "#238b45",'Exc L4-5 RORB GABRG1': "#41ab5d",'Exc L4-5 RORB IL1RAPL2': "#74c476",'Exc L5 ET': "#a1d99b",
        'Exc L5/6 IT Car3': "#c7e9c0",'Exc L5/6 NP': "#2ca25f",'Exc L5-6 RORB LINC02196': "#66c2a4",'Exc L6 CT': "#99d8c9",'Exc L6 THEMIS NFIA': "#ccece6",'Exc L6b': "#5aae61",'Exc NRGN': "#1b7837",'Exc RELN CHD7': "#d9f0d3",
        #Inhibitory neurons (blue)
        'Inh ALCAM TRPM3': "#08306b",'Inh CUX2 MSR1': "#08519c",'Inh ENOX2 SPHKAP': "#2171b5",'Inh FBN2 EPB41L4A': "#4292c6",'Inh GPC5 RIT2': "#6baed6", 'Inh L1 PAX6 CA4': "#9ecae1",
        'Inh L1-2 PAX6 SCGN': "#c6dbef",'Inh L1-6 LAMP5 CA13': "#9fb3c8",'Inh L3-5 SST MAFB': "#252b6c",'Inh L5-6 PVALB STON2': "#2c7fb8",'Inh L5-6 SST TH': "#3690c0",
        'Inh L6 SST NPY': "#74a9cf",'Inh LAMP5 NRG1 (Rosehip)': "#a6bddb",'Inh LAMP5 RELN': "#d0d1e6",'Inh PTPRK FAM19A1': "#1d3557",'Inh PVALB CA8 (Chandelier)': "#0570b0",
        'Inh PVALB HTR4': "#3690c0",'Inh PVALB SULF1': "#74a9cf",'Inh RYR3 TSHZ2': "#a6bddb",'Inh SGCD PDE3A': "#045a8d",'Inh SORCS1 TTN': "#2b8cbe",'Inh VIP ABI3BP': "#627d98",
        'Inh VIP CLSTN2': "#005082",'Inh VIP THSD7B': "#506680",'Inh VIP TSHZ2': "#41526a",
        #Astrocytes (orange)
        'Ast DPP10': "#d94801",'Ast GRM3': "#ff7f2a",'Ast TPST1': "#ffae66",
        #Microglia (red)
        'Mic DUSP1': "#660000",'Mic MKI67': "#990000",'Mic P2RY12': "#cc0000",'Mic TPT1': "#e60000",'T cells': "#ff1a1a",
        #Oligodendrocytes (pink)
        'Oli OPALIN': "#fa91c9",'Oli RASGRF1': "#ff4da6",
        #OPCs + COPs (brown)
        'OPC DOCK5': "#5c4033",'Opc GRIA4': "#8b5a2b",'Opc TPST1': "#d2b48c",
        #Vascular Cells (purple)
        'CAMs': "#3f007d",'Fib1': "#54278f",'Fib2': "#6a51a3",'Fib3': "#7b3294",'Per1': "#807dba",'Per2': "#9e9ac8",'aEndo': "#9e4db3",'aSMC': "#b084c1",'capEndo': "#cbc9e2",'vEndo': "#dadaeb",'vSMC': "#decbe4",
        #NA (black)
        'NA': "#000000"
        }
    condition_colors={"AD": "#8C6BB1", "ASD": "#1f77b4", "BD": "#41AB5D", "Control": "#E31A1C", "SZ": "#969696", "MS": "#FDB462"}
    datasource_colors={"Mathys": "#8C6BB1", "Wamsley": "#1f77b4", "Han": "#41AB5D", "Ruzicka": "#969696", "Macnair": "#FDB462"}
    if celltype_class_col:
        if not adata.obs[celltype_class_col].dtype.name=="category":
            adata.obs[celltype_class_col]=adata.obs[celltype_class_col].astype("category")
    if celltype_subtype_col:
        if not adata.obs[celltype_subtype_col].dtype.name=="category":
            adata.obs[celltype_subtype_col]=adata.obs[celltype_subtype_col].astype("category")
    if datasource_col:
        if not adata.obs[datasource_col].dtype.name=="category":
            adata.obs[datasource_col]=adata.obs[datasource_col].astype("category")
    if condition_col:
        if not adata.obs[condition_col].dtype.name=="category":
            adata.obs[condition_col]=adata.obs[condition_col].astype("category")
    if adata is None:
        if broad==False:
            return(subtype_colors)
        else:
            return(broad_colors)
    if celltype_class_col:
        colors=[broad_colors[cat] for cat in adata.obs[celltype_class_col].cat.categories]
        uns_classname=f"{celltype_class_col}_colors"
        adata.uns[uns_classname]=colors
    if celltype_subtype_col:
        colors=[subtype_colors[cat] for cat in adata.obs[celltype_subtype_col].cat.categories]
        uns_subname=f"{celltype_subtype_col}_colors"
        adata.uns[uns_subname]=colors
    return(adata)

def reorder_subtypes(celltypes):
    '''
    Used to define subtypes in priority order for graphing
    '''
    groups=[['exc'],['inh'],['ast'],['opc'], ['oli'],['mic', 'mg', 't cells'],['vas', 'endo', 'peri']]
    
    def get_priority(item):
        item_lower=item.lower()
        for i, keywords in enumerate(groups):
            if any(k in item_lower for k in keywords):
                return i
        return len(groups)

    sorted_list=sorted(celltypes, key=lambda x: (get_priority(x), x.lower()))
    return sorted_list

def basic_preprocess(adata, batch_col=None, scrub=None, human=True):
    '''
    Used to perform basic preprocessing of singlecell data from h5ad format.
    '''
    print('pre-filtering steps shape: ', adata.shape, "\n")
    if human==True:
        adata.var["mt"]=adata.var_names.str.startswith("MT-") #human mitochondrial genes specified
        adata.var["ribo"]=adata.var_names.str.startswith(("RPS", "RPL"))
        adata.var["hb"]=adata.var_names.str.contains("^HB[^(P)]")
    else:
        adata.var["mt"]=adata.var_names.str.startswith("mt-") #human mitochondrial genes specified
        adata.var["ribo"]=adata.var_names.str.startswith(("Rps", "Rpl"))
        adata.var["hb"]=adata.var_names.str.contains("^Hb[ab]")
        
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)
    print('number of cells with greater than 10% mitochondrial genes',len(adata[(adata.obs.pct_counts_mt > 10)]), "\n")
    adata=adata[(adata.obs.pct_counts_mt < 10)] #filtering out cells with greater than 10% mitochondrial genes
    #Ensuring that batch is category type
    if batch_col:
        if adata.obs[batch_col].dtype.name != 'category':
            adata.obs[batch_col]=adata.obs[batch_col].apply(lambda x: f'Batch{x}' if adata.obs[batch_col].dtype.name != 'category' else x).astype('category')
    #scrublet...
    if scrub:
        print('starting scrublet analysis')
        print('pre-scrublet adata shape: ', adata.shape, "\n")
        if batch_col:
            try:
                sc.pp.scrublet(adata, batch_key=batch_col)
                predicted_doublet_idx=list(adata[adata.obs.predicted_doublet==True].obs.index)
                print(len(predicted_doublet_idx), " number of doublets predicted across the dataset ","\n")
                if len(predicted_doublet_idx) > 0:
                    adata=adata[~adata.obs.index.isin(predicted_doublet_idx)] #filtering cells in predicted doublet index list out of adata object
                    print('doublets removed from', 'new shape: ', adata.shape, "\n")
                else:
                    print('no doublets to remove from', "\n")
            except Exception as e:
                print('scrublet with batch but without batch separation failed for ', "\n")
                print(f'Error: {str(e)}', "\n")
                batch_ids=adata.obs[batch_col].unique()
                predicted_doublet_idx=[]
                for batch in batch_ids:
                    print("scrubbing batch: ", batch)
                    adata_batch=adata[adata.obs["batch"]==batch]
                    sc.pp.scrublet(adata_batch) #performing scrublet on individual batches #may need to look at adjusting number of principal components later
                    predicted_doublet_idx.append(list(adata_batch[adata_batch.obs.predicted_doublet==True].obs.index)) #pulling a list of cell where predicted_doublet is true and adding this to the list of predicted doublet cells for the entire batch
    
                flat_predicted_doublet_idx=[cell_id for batch in predicted_doublet_idx for cell_id in batch]
                print(len(flat_predicted_doublet_idx), ": number of doublets predicted across the dataset ", "\n")     
                if len(flat_predicted_doublet_idx) > 0:
                    print('Done with scrublet, starting filter out of doublets \n')
                    adata=adata[~adata.obs.index.isin(flat_predicted_doublet_idx)] #filtering cells in predicted doublet index list out of adata object
                    print('done filtering of doublets from adata ', "\n")
                else: #no cells to filter out so continue
                    print('No predicted doublets. Adata shape: ', adata.shape, "\n")
        else: #batch being None
            sc.pp.scrublet(adata)
            predicted_doublet_idx=list(adata[adata.obs.predicted_doublet==True].obs.index)
            print(len(predicted_doublet_idx), " number of doublets predicted across the dataset ", '\n')
            adata=adata[~adata.obs.index.isin(predicted_doublet_idx)] #filtering cells in predicted doublet index list out of adata object
            print('done filtering of doublets from adata ')
        print('adata shape after scrublet: ', adata.shape, '\n')
    #normalize data...
    if adata.X.max()<20:
        print('Already logged and normalized counts \n')
    else:
        print(adata, 'normalizing and log transforming... \n')
        adata.layers["counts"]=adata.X.copy()
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    #HV genes
    if batch_col:
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=batch_col)
    else:
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    #clustering
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    return(adata)
