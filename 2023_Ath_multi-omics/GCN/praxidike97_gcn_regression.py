<<<<<<< HEAD
""" 
GCN for Classification from praxidike97 on GitHub adapted for regression with genotype data.

GCN source: https://github.com/praxidike97/GraphNeuralNet/blob/master/main.py
layer source: https://github.com/dizhu-gis/SRGCNN/blob/main/SRGCNN_demo.ipynb
Adjacency matrix source: https://github.com/txWang/MOGONET

To run on command line if CUDA error (device-side assert triggered): 
    CUDA_LAUNCH_BLOCKING=1 python praxidike97_gcn_regression.py 

Data used:
    Planetoid: 
        Data(x=[2708, 1433], edge_index=[2, 10556], y=[2708], train_mask=[2708], val_mask=[2708], test_mask=[2708])
        num_node_features = 1433; num_classes = 7
    
    Genotype: Data(x=[383, 1000], edge_index=[2, 5744], y=[383], train_mask=[383], test_mask=[383])
        1000 features, 383 instances
"""

=======
# GCN for Classification from praxidike97 on GitHub adapted for regression
# Source: https://github.com/praxidike97/GraphNeuralNet/blob/master/main.py
# To run on command line if error: CUDA_LAUNCH_BLOCKING=1 python multi-omics/GCN/praxidike97_gcn.py 
""" 
Classification: run on Planetoid data
Info: 
- Data(x=[2708, 1433], edge_index=[2, 10556], y=[2708], train_mask=[2708], val_mask=[2708], test_mask=[2708])
- num_node_features = 1433; num_classes = 7

Regression: run on genotype data
"""
from torch_geometric.data.in_memory_dataset import InMemoryDataset
from torch_geometric.datasets import Planetoid
>>>>>>> 314110fd1adc175794a25433dd6fdad0398f8783
import torch
from torch import nn
from torch import tensor
import torch.nn.functional as F
from torch_geometric.datasets import Planetoid
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops, degree
<<<<<<< HEAD
from torch_geometric.data import Data
=======
from torch_geometric.data import InMemoryDataset, Data
#from torch_geometric.nn import GCNConv # graph convolutional layer
>>>>>>> 314110fd1adc175794a25433dd6fdad0398f8783
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import datatable as dt

<<<<<<< HEAD
# Graph Convolutional layer
=======
# Custom graph convolutional layer
>>>>>>> 314110fd1adc175794a25433dd6fdad0398f8783
class GCNConv(MessagePassing):
    def __init__(self, in_channels, out_channels):
        super(GCNConv, self).__init__(aggr='add')  # "Add" aggregation
        self.lin = torch.nn.Linear(in_channels, out_channels)

    def forward(self, x, edge_index):
        # Step 1: Add self-loops
        edge_index, _ = add_self_loops(edge_index, num_nodes=x.size(0))

        # Step 2: Multiply with weights
        x = self.lin(x)

        # Step 3: Calculate the normalization
        row, col = edge_index
<<<<<<< HEAD
        deg = degree(row, x.size(0), dtype=x.dtype)
=======
        deg = degree(row, x.size(0), dtype=x.dtype) ####### this step here is not working with genotype data, but everything works for planetoid data
>>>>>>> 314110fd1adc175794a25433dd6fdad0398f8783
        deg_inv_sqrt = deg.pow(-0.5)
        norm = deg_inv_sqrt[row] * deg_inv_sqrt[col]

        # Step 4: Propagate the embeddings to the next layer
        return self.propagate(edge_index, size=(x.size(0), x.size(0)), x=x,
                              norm=norm)

    def message(self, x_j, norm):
        # Normalize node features.
        return norm.view(-1, 1) * x_j

<<<<<<< HEAD
# Graph convolutional layer for regression
class GCNConv_reg(nn.Module):
    """Basic graph convolution operation that incorporate both spatial lagged X 
    and spatial lagged Y (to be used in the basic SRGCNNs model)"""
    def __init__(self, f_in, f_out, use_bias=True, activation=F.relu):#hidden layer with relu activation
        super().__init__()
        self.f_in = f_in
        self.f_out = f_out
        self.use_bias = use_bias
        self.activation = activation
        self.weight = nn.Parameter(torch.FloatTensor(f_in, f_out),requires_grad=True)###requires_grad:whether allow weights to be updated
        self.bias = nn.Parameter(torch.FloatTensor(f_out)) if use_bias else None
        self.initialize_weights()
    
    def initialize_weights(self):
        nn.init.constant_(self.weight,1)
        if self.use_bias: nn.init.constant_(self.bias,0)
        
    def forward(self, input, adj):
        support = torch.mm(input, self.weight)
        output = torch.mm(adj, support) #adj here has to be renormalized     
        
        if self.use_bias: output.add_(self.bias)
        if self.activation is not None: output=self.activation(output) 
        
        return output

# GCN Model
=======
# GCN Model two-layers
>>>>>>> 314110fd1adc175794a25433dd6fdad0398f8783
class Net(torch.nn.Module):
    def __init__(self, dataset):
        super(Net, self).__init__()
        # for classification
<<<<<<< HEAD
        self.conv1 = GCNConv(dataset.num_node_features, 16) 
        self.conv2 = GCNConv(16, dataset.num_classes)
        # for regression
        #self.conv1 = GCNConv(1000, 400) 
        #self.conv2 = GCNConv(400, 383)
=======
        #self.conv1 = GCNConv(dataset.num_node_features, 16)
        #self.conv2 = GCNConv(16, dataset.num_classes)
        # for regression
        self.conv1 = GCNConv(1000, 400)
        self.conv2 = GCNConv(400, 319)
>>>>>>> 314110fd1adc175794a25433dd6fdad0398f8783

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        #x = x.to(device)
        #edge_index = edge_index.type(torch.cuda.FloatTensor).cuda(device)
        #print('x', x, x.size())
        #print('edge_index', edge_index, edge_index.size())
        x = self.conv1(x, edge_index)
        x = F.relu(x) # this is messing up the structure of the tensor for genotype data!! :( but it's fine for planetoid. Why is that?
        x = F.dropout(x, training=self.training)
        #print('dropout x', x, x.size())
        x = self.conv2(x, edge_index)
<<<<<<< HEAD
        #print('x 2', x, x.size())
        
        #return x # for regression
        return F.log_softmax(x, dim=1) # for classification
        
=======

        return x #F.log_softmax(x, dim=1) # softmax paired with nll_loss function for classification
>>>>>>> 314110fd1adc175794a25433dd6fdad0398f8783

# Graph of dataset
def plot_dataset(dataset):
    edges_raw = dataset.data.edge_index.numpy()
    edges = [(x, y) for x, y in zip(edges_raw[0, :], edges_raw[1, :])]
    labels = dataset.data.y.numpy()

    G = nx.Graph()
    G.add_nodes_from(list(range(np.max(edges_raw))))
    G.add_edges_from(edges)
    plt.subplot(111)
    options = {
                'node_size': 30,
                'width': 0.2,
    }
    nx.draw(G, with_labels=False, node_color=labels.tolist(), cmap=plt.cm.tab10, font_weight='bold', **options)
    plt.savefig("planetoid_graph.png")


# Evaluate the model on the test set
def test(data, train=True):
    model.eval()

    correct = 0
    pred = model(data).max(dim=1)[1]
    #print("TEST-------------")
    #print('pred', pred, pred.size())
    #print('data.train_mask', data.train_mask, data.train_mask.size())
    #print('data.y', data.y, data.y.size())
    #sprint('data.y[data.train_mask]', data.y[data.train_mask])
    if train:
        correct += pred[data.train_mask].eq(data.y[data.train_mask]).sum().item()
        return correct / (len(data.y[data.train_mask]))
    else:
        correct += pred[data.test_mask].eq(data.y[data.test_mask]).sum().item()
        return correct / (len(data.y[data.test_mask]))

# train the model
def train(data, plot=False):    
    train_accuracies, test_accuracies = list(), list()
    start = time.time()
    for epoch in range(100):
        model.train()
        optimizer.zero_grad()
        #print('data', data, data.size())
        out = model(data)
        # IndexError: The shape of the mask [383] at index 0 does not match the shape of the indexed tensor [2, 7] at index 0
        #print('train mask', data.train_mask, data.train_mask.size())
        #print('out', out, out.size())
        #print('out[data.train_mask]', out[data.train_mask], out[data.train_mask].size())
        #print('data.y', data.y, data.y.size())
        #print('data.y[data.train_mask]', data.y[data.train_mask], data.y[data.train_mask].size())
        
        # for regression
        #loss = torch.nn.MSELoss()
        #output = loss(out[data.train_mask], data.y[data.train_mask]) 
        #output.backward()

        # for classification
        loss = F.nll_loss(out[data.train_mask], data.y[data.train_mask])
        loss.backward()
        optimizer.step()


        train_acc = test(data)
        test_acc = test(data, train=False)

        train_accuracies.append(train_acc)
        test_accuracies.append(test_acc)
        print('Epoch: {:03d}, Loss: {:.5f}, Train Acc: {:.5f}, Test Acc: {:.5f}'.
              format(epoch, loss, train_acc, test_acc))
    end = time.time()
    print("Elapsed time: ", end-start)

    if plot: # plot AUC curve
        plt.plot(train_accuracies, label="Train accuracy")
        plt.plot(test_accuracies, label="Validation accuracy")
        plt.xlabel("# Epoch")
        plt.ylabel("Accuracy")
        plt.title("Planetoid Dataset")
        plt.legend(loc='upper right')
        plt.savefig("auc_praxidike_gcn.png")

'''
MOGONET Method for Computing Adjacency Matrix
1. Calculate adjacency matrix parameter (cal_adj_mat_parameter)
    a. Compute cosine similarity matrix (cosine_distance_torch)
2. Generate adjacency matrix tensor (gen_adj_mat_tensor)
    a. Compute cosine similarity matrix (cosine_distance_torch)
    b. Generate graph from distance matrix (graph_from_dist_torch)
    c. Compute identity matrix
    d. Compute adjacency matrix (not equation 2 in paper)
3. Generate test adjacency matrix tensor (gen_test_adj_mat_tensor)
'''
def cosine_distance_torch(x1, x2=None, eps=1e-8):
    x2 = x1 if x2 is None else x2
    w1 = x1.norm(p=2, dim=1, keepdim=True)
    w2 = w1 if x2 is x1 else x2.norm(p=2, dim=1, keepdim=True)
    return 1 - torch.mm(x1, x2.t()) / (w1 * w2.t()).clamp(min=eps)

def cal_adj_mat_parameter(edge_per_node, data, metric="cosine"):
    assert metric == "cosine", "Only cosine distance implemented"
    dist = cosine_distance_torch(data, data)
    parameter = torch.sort(dist.reshape(-1,)).values[edge_per_node*data.shape[0]]
    return parameter.data.cpu().numpy().item()

def graph_from_dist_tensor(dist, parameter, self_dist=True):
    if self_dist:
        assert dist.shape[0]==dist.shape[1], "Input is not pairwise dist matrix"
    g = (dist <= parameter).float() # binary 0s and 1s for False/True
    if self_dist:
        diag_idx = np.diag_indices(g.shape[0])
        g[diag_idx[0], diag_idx[1]] = 0
    return g

def to_sparse(x):
    x_typename = torch.typename(x).split('.')[-1]
    sparse_tensortype = getattr(torch.sparse, x_typename)
    indices = torch.nonzero(x)
    if len(indices.shape) == 0:  # if all elements are zeros
        return sparse_tensortype(*x.shape)
    indices = indices.t()
    values = x[tuple(indices[i] for i in range(indices.shape[0]))]
    return sparse_tensortype(indices, values, x.size()), indices, values

def gen_adj_mat_tensor(data, parameter, metric="cosine"):
    assert metric == "cosine", "Only cosine distance implemented"
    dist = cosine_distance_torch(data, data)
    g = graph_from_dist_tensor(dist, parameter, self_dist=True)
    if metric == "cosine":
        adj = 1-dist
    else:
        raise NotImplementedError
    adj = adj*g 
    adj_T = adj.transpose(0,1)
    I = torch.eye(adj.shape[0])
    cuda = True if torch.cuda.is_available() else False
    adj = adj + adj_T*(adj_T > adj).float() - adj*(adj_T > adj).float()
    adj = F.normalize(adj + I, p=1)
    adj = to_sparse(adj)   
    return adj

# create PyTorch Data object using genotype data
def test_geno():
    geno_data="SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv"
    pheno_data="Phenotype_value_383_common_accessions_2017_Grimm.csv"
    test_mask="test_20perc.txt"
    test_mask = pd.read_csv(test_mask, header=None)
    geno = dt.fread(geno_data) # read in genotype data
    geno = geno.to_pandas() # convert dataframe to pandas dataframe
    geno = geno.sort_values(by=geno.columns[0], axis=0) # sort values by sample ID
    geno = geno.set_index(geno.columns[0], drop=True) # set index to sample ID
    geno_sub = geno.iloc[:,0:1000] # 383 samples rows
    features = geno_sub.columns # columns as features
    
    pheno = pd.read_csv(pheno_data, index_col=0) # read in phenotype data
    label = pheno.FT10_mean
    
    
    # Split geno and pheno into training and testing sets
    X_train = geno_sub.loc[~geno_sub.index.isin(test_mask[0])]
    X_test = geno_sub.loc[geno_sub.index.intersection(test_mask[0])]
    y_train = pheno.loc[~pheno.index.isin(test_mask[0])]
    y_test = pheno.loc[pheno.index.intersection(test_mask[0])]
    #y_train = np.random.randint(0,6,306) # attempting classification
    #y_test = np.random.randint(0,6,77)
    #y = np.random.randint(0,6,383)

    # Create masks
    train_mask = tensor([i in np.array(X_train.index) for i in np.array(geno.index)])
    test_mask = tensor([i in np.array(X_test.index) for i in np.array(geno.index)])

    # Convert to PyTorch tensors
    geno_sub = torch.tensor(geno_sub.values.astype(np.float32))
    geno_sub[geno_sub==-1] = 0 # convert -1 to 0, just in case
    label = torch.tensor(label.values.astype(np.float32))
    X_train = torch.tensor(X_train.values.astype(np.float32))
    X_test = torch.tensor(X_test.values.astype(np.float32))
    y_train = torch.tensor(y_train)
    y_test = torch.tensor(y_test)
    #y = torch.tensor(y)

    # Compute adjacency matrix
    adj_parameter = 2 # edge_per_node
    adj_parameter_adaptive = cal_adj_mat_parameter(adj_parameter, X_train.t(), "cosine")
    adj_train, edge_index, adj_values = gen_adj_mat_tensor(X_train.t(), adj_parameter_adaptive, "cosine")

    # Create PyTorch geometric Data object
    return Data(x=geno_sub, edge_index=edge_index, y=label, train_mask=train_mask, test_mask=test_mask)

if __name__ == "__main__":
    import gc
    gc.collect()
    torch.cuda.empty_cache()
    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu') # check device
    print(device)
    t = torch.cuda.get_device_properties(0).total_memory
    r = torch.cuda.memory_reserved(0)
    a = torch.cuda.memory_allocated(0)
    f = r-a  # free inside reserved
    print(t, r, a, f)

<<<<<<< HEAD
    # planetoid data for classification
    dataset = Planetoid(root='/tmp/Cora', name='Cora') # load dataset (2708 nodes/input tensor vars, 1433 instances)
    data = dataset[0].to(device)
    
    # genotype data for regression
    #dataset = test_geno()
    #data = dataset.cuda(device)
    #print("data", data)
=======
    # Classification: planetoid data
    #dataset = Planetoid(root='/tmp/Cora', name='Cora') # load dataset (2708 nodes/input tensor vars, 1433 instances)
    #data = dataset[0].to(device)
    
    # Regression: genotype data
    dataset = test_geno()
    data = dataset.cuda(device)
>>>>>>> 314110fd1adc175794a25433dd6fdad0398f8783

    model = Net(dataset) # create GCN model
    model.cuda(device) # send to gpu
    
    # Optimizer    
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)
    # Train the model
    train(data, plot=True)
