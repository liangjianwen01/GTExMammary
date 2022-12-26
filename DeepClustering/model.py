import torch
import torch.nn as nn
import torch.nn.functional as F

class MLPclassifyer(nn.Module):
    def __init__(self,in_feature,out_feature):
        super(MLPclassifyer,self).__init__()
        
        self.hidden1=nn.Linear(in_features=in_feature,out_features=100,bias=True)#in_feature*100 in_feature个属性特征
       
        self.hidden2=nn.Linear(100,100)#100*100
        
        self.hidden3=nn.Linear(100,50)#100*50
       
        self.predict=nn.Linear(50,out_feature)#50*1  
        
        self.norm=nn.LayerNorm(in_feature)
        self.norm2=nn.LayerNorm(50)
        self.softmax = nn.Softmax(dim=1)
    def forward(self,x):
        x=self.norm(x)
        x=F.relu(self.hidden1(x))
        x=F.relu(self.hidden2(x))
        x=F.relu(self.hidden3(x))
        x=self.norm2(x)
        output=self.softmax(self.predict(x))
        #output = self.predict(x)
        return output,x


