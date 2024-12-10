########################## file information #############################
# Author: Pengfei Ren
# Date: 2024-12-04 21:48:39
# LastEditTime: 2024-12-04 23:32:58
# Description: 
# FilePath: /undefined/Users/morsouron/Desktop/Project/benchmark/code/8_st_annotation.py
#########################################################################

# selina
sys.path.append('../../')
from Selina.selina import *
def predict(train_data,train_meta,adata):
    #common genes
    common_genes = np.intersect1d(train_data.var_names.tolist(),adata.var_names.tolist())
    print(len(common_genes))
    train_data = train_data[:,common_genes].to_df().T
    adata = adata[:,common_genes.tolist()]
    # train
    train_data, celltypes, platforms, genes = preprocessing(train_data,train_meta)
    print(Counter(celltypes))
    ct_dic, plat_dic = label2dic(celltypes), label2dic(platforms)
    nfeatures, nct, nplat = train_data.shape[0], len(ct_dic), len(plat_dic)
    network = train(train_data, params_train, celltypes, platforms, nfeatures, nct, nplat, ct_dic, plat_dic, device)
    # predict
    network = Autoencoder(network, nfeatures, nct)
    network = Normal_Classifier(network).to(device) 
    pred_labels, pred_prob = test(adata, network, ct_dic, device)
    return pred_labels

#spoint
def predict(sc_ad,st_ad):
    spoint_model = Spoint.init_model(sc_ad,st_ad,celltype_key='major',deg_method='t-test',sm_size=100000,use_gpu=True)
    spoint_model.train(max_steps=5000, batch_size=1024,save_mode='best')
    st_data = torch.tensor(spoint_model.st_data).to('cpu')
    model = spoint_model.model.to('cpu')
    model.eval()
    latent, pre, _ = model(st_data)
    pre = pre.detach().numpy()
    pre = pd.DataFrame(pre,columns=spoint_model.clusters,index=spoint_model.st_ad.obs_names)
    max_columns = pre.idxmax(axis=1)
    st_ad.obs['ct'] = 'Unknown'
    st_ad.obs.loc[max_columns.index,'ct'] = max_columns.values
    return(st_ad.obs['ct'].values)

#tangram
def tangram_pred(st_ad,sc_ad,tissue):
    df_genes = pd.read_csv('/home/renpf/benchmark/res/degs/scrna/' + tissue + '_top50.txt',header=0,sep='\t')
    markers = df_genes['gene'].tolist()
    tg.pp_adatas(sc_ad, st_ad, genes=markers)
    assert sc_ad.uns['training_genes'] == st_ad.uns['training_genes']
    ad_map = tg.map_cells_to_space(
        adata_sc=sc_ad,
        adata_sp=st_ad,
        device='cuda:1',
        mode = "clusters",
        cluster_label='major'
    )
    tg.project_cell_annotations(ad_map, st_ad, annotation='major')
    return(st_ad.obsm['tangram_ct_pred'])

#tacco
for i in range(len(tissues)):
    tissue = tissues[i]
    sc_ad = ref[i]
    sc_ad.X = sc_ad.X.astype(np.float32)
    st_ad = sc.read_h5ad('st.h5ad')
    st_ad.X = st_ad.X.astype(np.float32)
    tc.tl.annotate(st_ad,sc_ad,'major',result_key='major')

# celltypist
def pred(sc_ad,st_ad,tissue):
    sc.pp.normalize_total(sc_ad, target_sum=1e4)
    sc.pp.log1p(sc_ad)
    sc.pp.normalize_total(st_ad, target_sum=1e4)
    sc.pp.log1p(st_ad)
    model = celltypist.train(sc_ad, labels = 'major', n_jobs = 10, feature_selection = True)
    predictions = celltypist.annotate(st_ad, model = model, majority_voting = False)
    return(predictions)

# merge annotations generated with different tools
def annotation(adata,directory_path,seg_type=None):
    files = os.listdir(directory_path)
    if seg_type is not None:
        files = [file for file in files if seg_type in file]
    data = []
    for file in files:
        file_path = os.path.join(directory_path, file)
        data.append(pd.read_csv(file_path, header=None))
    data = pd.concat(data,axis=1)
    data['final_annotation'] = data.apply(vote, axis=1)
    intersect = []
    for i in range(5):
        intersect.append(np.sum(data['final_annotation']==data.iloc[:,i]))
    method_index = intersect.index(max(intersect))
    data.loc[data['final_annotation']=='unknown','final_annotation'] = data.iloc[:,method_index][data['final_annotation']=='unknown']
    final_annotation = data['final_annotation'].tolist()
    return(final_annotation)