import os
import time
import numpy as np
import pandas as pd 
import datatable as dt
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch_geometric.datasets import QM7b
from torch_geometric.datasets import Planetoid
os.chdir("/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/multi-omics/GCN")
import simple_GCN

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

def cosine_distance_torch(x1, x2=None, eps=1e-8):
    x2 = x1 if x2 is None else x2
    w1 = x1.norm(p=2, dim=1, keepdim=True)
    w2 = w1 if x2 is x1 else x2.norm(p=2, dim=1, keepdim=True)
    return 1 - torch.mm(x1, x2.t()) / (w1 * w2.t()).clamp(min=eps)

def cal_adj_mat_parameter(edge_per_node, data, metric="cosine"):
    assert metric == "cosine", "Only cosine distance implemented"
    dist = cosine_distance_torch(data, data)
    parameter = torch.sort(dist.reshape(-1,)).values[edge_per_node*data.shape[0]]
    return np.asscalar(parameter.data.cpu().numpy())

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
    if cuda:
        I = I.cuda('cuda:0')
    adj = adj + adj_T*(adj_T > adj).float() - adj*(adj_T > adj).float()
    adj = F.normalize(adj + I, p=1)
    adj = to_sparse(adj)   
    return adj

# Try running Planetoid data in GCN
def test_gcn(dataset, data):
    # Define the network
    net = simple_GCN.GCN(in_channels = dataset.num_node_features, out_channels = [100, 16, dataset.num_classes])
    net = net.to(device='cuda:0')
    print(net)

    optimizer = torch.optim.SGD(net.parameters(), lr=0.2) # optimizer (variant of gradient descent)
    loss_func = torch.nn.MSELoss()  # mean squared loss function for regression

    start = time.time()
    # Train the network
    print("Training... ")
    for epoch in range(500):
        current_loss = 0.0
        print(f"\nTest: Epoch {epoch}")
        optimizer.zero_grad()                   # clear gradients for next train
        prediction = net(data.x[data.train_mask], data.edge_index) # prediction based on input x
        #loss = loss_func(prediction[data.train_mask], data.y[data.train_mask]) # for regression
        loss = F.nll_loss(prediction[data.train_mask], data.y[data.train_mask]) # for classification    # must be (1. nn output, 2. target)
        loss.backward()                         # backpropagation, compute gradients
        optimizer.step()                        # apply gradients

        # print statistics
        current_loss += loss.item()
        print('Loss after mini-batch %5d: %.3f' % (input + 1, current_loss / 500))
        current_loss = 0.0
    
    end = time.time()   
    elapsed = end - start
    print("Training complete; Elapsed time: ", elapsed)

    #dataset = QM7b(root="/tmp/Cora")
    #edges_raw = dataset.data.edge_index.numpy()
    #edges = [(x, y) for x, y in zip(edges_raw[0, :], edges_raw[1, :])]
    #labels = dataset.data.y.numpy()

def main():
    torch.cuda.empty_cache()
    cuda = True if torch.cuda.is_available() else False
    
    # Run Planetoid data through GCN
    dataset = Planetoid(root='/tmp/Cora', name='Cora') # load dataset (2708 nodes/input tensor vars, 1433 instances)
    data = dataset[0].to(device='cuda:0')
    test_gcn(dataset, data) # test dataset

    '''
    """ Run on genotype data """
    geno_sub, X_train, y_train, train_loader, test_loader, test_tensor = simple_GCN.split_data()

    # compute adjacency matrix
    adj_parameter = 2 # edge_per_node
    adj_parameter_adaptive = cal_adj_mat_parameter(adj_parameter, X_train.t(), "cosine")
    adj_train, adj_indices, values = gen_adj_mat_tensor(X_train.t(), adj_parameter_adaptive, "cosine")
    
    # Define the network
    net = simple_GCN.GCN(in_channels = X_train.size()[1], out_channels = [800, 500, X_train.size()[0]])
    net = net.to(device="cuda:0")
    print(net)

    optimizer = torch.optim.SGD(net.parameters(), lr=0.2) # optimizer (variant of gradient descent)
    loss_func = torch.nn.MSELoss()  # mean squared loss function for regression

    start = time.time()
    # Train the network
    print("Training... ")
    for epoch in range(500):
        current_loss = 0.0
        print(f"\nTest: Epoch {epoch}")
        #for input, batch in enumerate(train_loader, 0): # Data loader
        optimizer.zero_grad()                   # clear gradients for next train
        #inputs, target = batch
        print(X_train.dtype, adj_indices.dtype)
        prediction = net(X_train, adj_indices)#inputs)                # prediction based on input x
        print("prediction", prediction)
        print("y_train", y_train, y_train.size())
        loss = loss_func(prediction, y_train)#, target)    # must be (1. nn output, 2. target)
        print("loss", loss)
        loss.backward()                         # backpropagation, compute gradients
        optimizer.step()                        # apply gradients

        # print statistics
        current_loss += loss.item()
        print('Loss after mini-batch %5d: %.3f' % (input + 1, current_loss / 500))
        current_loss = 0.0
    
    end = time.time()   
    elapsed = end - start
    print("Training complete; Elapsed time: ", elapsed)

    # Evaluate trained model on test set
'''   

if __name__ == "__main__":
    main()