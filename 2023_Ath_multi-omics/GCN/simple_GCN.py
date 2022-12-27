import os
import time
import numpy as np
import pandas as pd 
import datatable as dt
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.metrics.pairwise import cosine_similarity
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils.data as data
import torch.optim as optim
from torch_geometric.nn import GCNConv
import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor, ModelCheckpoint
from torch_geometric.nn.conv.message_passing import MessagePassing

# Set seed for reproducibility
#pl.seed_everything(42)

# Ensure that all operations are deterministic on GPU for reproducibility
torch.backends.cudnn.determinstic = True
torch.backends.cudnn.benchmark = False

print(torch.version.cuda)
device = torch.device("cuda:0") if torch.cuda.is_available() else torch.device("cpu") # set to gpu or cpu
print(device)
torch.cuda.device_count() # number of gpus

if device.type == 'cuda':
    print(torch.cuda.get_device_name(0))
    print('Memory Usage:')
    print('Allocated:', round(torch.cuda.memory_allocated(0)/1024**3,1), 'GB')
    print('Cached:   ', round(torch.cuda.memory_cached(0)/1024**3,1), 'GB')

# Set path to directory containing data
PATH = "/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic"
os.chdir(PATH)


class GCN_layer(nn.Module): # source: MOGONET
    """
    A Graph Convolution Layer (GCN).
    """
    def __init__(self, in_channels, out_channels, bias: bool = True):
        super(GCN_layer, self).__init__()
        self.in_channels = in_channels # number of input nodes
        self.out_channels = out_channels # number of output nodes
        self.weight = nn.Parameter(torch.FloatTensor(in_channels, out_channels)) # initialize weights
        if bias:
            self.bias = nn.Parameter(torch.FloatTensor(out_channels))
        nn.init.xavier_normal_(self.weight.data)
        if self.bias is not None:
            self.bias.data.fill_(0.0)

    def forward(self, x, adj): #data
        support = torch.mm(x, self.weight)
        print("support", support, support.size(), support.dtype)
        print("adj", adj, adj.size(), adj.dtype)
        output = torch.sparse.mm(adj.float(), support.float())
        if self.bias is not None:
            return output + self.bias
        else:
            return output

class GCN(nn.Module):
    """
    A three-layer GCN.
    """
    def __init__(self, in_channels, out_channels): # size of input, size of outputs
        super(GCN, self).__init__()
        torch.manual_seed(42) # for reproducibility
        # Add layers
        self.layer1 = GCN_layer(in_channels, out_channels[0], bias = True)#, cached = True) #GCNConv
        self.layer2 = GCN_layer(out_channels[0], out_channels[1], bias = True)#, cached = True) #GCNConv
        self.layer3 = GCN_layer(out_channels[1], out_channels[2], bias = True)#, cached = True) #GCNConv
        self.classifier = nn.Sequential(nn.Linear(in_channels, out_channels[2]))

    def forward(self, x, adj):
        print("x", x, x.size(), x.dtype)
        print("adj", adj, adj.size(), adj.dtype)
        x = self.layer1(x, adj) # graph convolutional layer
        print("layer1", x, x.size())
        x = F.leaky_relu(x, 0.25) # activation function
        print("leaky_relu",x, x.size())
        x = F.dropout(x, training = self.training) # regularization
        print("regularization", x, x.size())
        x = self.layer2(x, adj)
        print("layer2", x, x.size())
        x = F.leaky_relu(x, 0.25)
        print(x)
        x = F.dropout(x, training=self.training)
        print(x)
        x = self.layer3(x, adj)
        print("layer3", x, x.size())
        x = F.leaky_relu(x, 0.25)
        print("leaky_relu:", x, x.size())
        # Apply the linear classifier
        #x = self.classifier(x)
        #x = F.log_softmax(x, dim=1)
        print("out: ", x, x.size())
        return x
 
 # include function to split data

def split_data(): # geno_data, pheno_data args; move to simple_GCN
    geno_data="SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv"
    pheno_data="Phenotype_value_383_common_accessions_2017_Grimm.csv"

    geno = dt.fread(geno_data) # read in genotype data
    geno = geno.to_pandas() # convert dataframe to pandas dataframe
    geno = geno.sort_values(by=geno.columns[0], axis=0) # sort values by sample ID
    geno = geno.set_index(geno.columns[0], drop=True) # set index to sample ID
    geno_sub = geno.iloc[:,0:1000]
    features = geno_sub.columns # columns as features

    pheno = pd.read_csv(pheno_data, index_col=0) # read in phenotype data
    label = pheno.FT10_mean # flowering time as label

    # Split geno and pheno into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(geno_sub, label, test_size=64)

    # Convert to PyTorch tensors
    geno_sub = torch.tensor(geno_sub.values.astype(np.float32))
    label = torch.tensor(label.values.astype(np.float32))
    X_train = torch.tensor(X_train.values.astype(np.float32), device = "cuda:0")
    X_test = torch.tensor(X_test.values.astype(np.float32), device = "cuda:0")
    y_train = torch.tensor(y_train.values.astype(np.float32), device = "cuda:0")
    y_test = torch.tensor(y_test.values.astype(np.float32), device = "cuda:0")
    
    train_tensor = data.TensorDataset(X_train, y_train) # Tensor training dataset to feed; add to cuda device in loader?
    train_loader = data.DataLoader(dataset = train_tensor, batch_size = 1, shuffle = True, pin_memory = False) # data loader to feed training data to network
    test_tensor = data.TensorDataset(X_test, y_test)
    test_loader = data.DataLoader(dataset=test_tensor, shuffle = True, pin_memory = False)
    return geno_sub, X_train, y_train, train_loader, test_loader, test_tensor

# Thought: Do I need to normalize my data?
# Does the input channel need to be a row vector or a matrix?