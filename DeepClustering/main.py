import torch
import Cluster
from sklearn.preprocessing import StandardScaler
import pandas as pd
import torch.utils.data as Data
import numpy as np
from sklearn import preprocessing

def loadcsv(file):
    source_data = pd.read_csv(file,index_col=0)
    data_list = source_data.values.tolist()
    return data_list

def KMM(data,n_cluster):
    kmm = Cluster.KMeans(n_clusters=n_cluster).fit(data)
    label = kmm.labels_
    return label

def data_process(data,label,batch_size):
    # scale=StandardScaler()
    # X_train = scale.fit_transform(data)
    train_xt=torch.from_numpy(data.astype(np.float32))
    train_yt=torch.from_numpy(label.astype(np.compat.long))
    train_data = Data.TensorDataset(train_xt,train_yt)
    train_loader = Data.DataLoader(dataset=train_data,batch_size=batch_size,shuffle=True,num_workers=0)
    return train_loader

if __name__ == '__main__':
    lose_diary = {}
    parameter_diary = {}
    nclass = 4 
    source_data = loadcsv(file='feature.csv')
    N_data = len(source_data)
    source_data = np.array(source_data)

    label_Y = KMM(preprocessing.scale(source_data),nclass)#改类别数量
    label_Y = np.array(label_Y)

    data_for_MLP = data_process(data=source_data,label=label_Y,batch_size=12)

    mlp_data = torch.from_numpy(source_data.astype(np.float32))

    cluster_result = []
    cluster_result.append(label_Y)

    for i in range(50):
        print(i)
        model = Cluster.creat_MLP_model(mlp_data[0].shape[0],nclass)
        parameter_diary['epoch'+str(i)],lose_diary['epoch'+str(i)]=model.train(train_loader=data_for_MLP,epo=500)
        label_Y,new_data = model.MLP_Cluster(data=mlp_data,ncluster=nclass)
        mlp_data = new_data
        data_for_MLP = data_process(data=np.array(mlp_data.tolist()),label=label_Y,batch_size=16)
        cluster_result.append(label_Y)
    cluster_df = pd.DataFrame({'id':list(range(1,N_data+1,1)),'cluster':cluster_result[-1]})
    cluster_df.to_csv('cluster_result.csv')
    






