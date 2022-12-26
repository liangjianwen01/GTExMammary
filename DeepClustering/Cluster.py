import model
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from torch.optim import Adam



class creat_MLP_model():
    def __init__(self,n,n2):
        self.mlp = model.MLPclassifyer(in_feature=n,out_feature=n2)
        self.optimizer = Adam(self.mlp.parameters(),lr=0.001,betas=(0.9,0.999),eps=1e-08,weight_decay=0,amsgrad=False)
        self.loss = model.nn.CrossEntropyLoss()
        self.train_loss_all=[]
        # self.index = 0
        self.parameters = []
        self.cluster = []
        

    def train(self,train_loader,epo):
        for epoch in range(epo):
            train_loss=0
            train_num=0
            for step,(b_x,b_y) in enumerate(train_loader):
                output,_=self.mlp.forward(b_x)
                loss=self.loss(output,b_y)
                self.optimizer.zero_grad()
                loss.backward()
                self.optimizer.step()
                train_loss+=loss.item()*b_x.size(0)
                train_num+=b_x.size(0)
            self.train_loss_all.append(train_loss/train_num)
            # self.index += 1
            # if self.index==10:
            #     self.index = 0
            if (epoch>0 and (self.train_loss_all[epoch] < self.train_loss_all[epoch-1])):
                self.parameters.append(self.mlp.parameters())
        return self.parameters,self.train_loss_all

    def MLP_Cluster(self,data,ncluster):
        _,new_feature = self.mlp.forward(data)
        kmm = KMeans(n_clusters=ncluster).fit(new_feature.tolist())
        self.cluster=kmm.labels_
        self.centroid = kmm.cluster_centers_
        return self.cluster,new_feature
            





