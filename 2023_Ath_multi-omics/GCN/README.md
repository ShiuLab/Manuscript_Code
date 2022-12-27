<<<<<<< HEAD
# File Descriptions
|Files|Description|
|-----|-----------|
|calc_adj.py|Calculate adjacency matrix|
|main.py|Main program to run simple_GCN.py|
|praxidike97_gcn_regression.py|Adapted GCN for classification from [praxidike97](https://github.com/praxidike97/GraphNeuralNet/blob/master/main.py) to do regression on genotype data|
|pycache|Compiled simple_GCN.py to import in main.py|
|simple_GCN.py|GCN based on [MOGONET](https://github.com/txWang/MOGONET) for regression on genotype data|
|SNP_Dataset.py|Attempt to make a PyTorch geometric dataset using genotype data|
=======
Files in this directory:
------------------------
- __pycache__ : contains the compiled simple GCN code
- __main.py__ : implements simple GCN on genotype data
- __praxidike97_gcn_regression.py__ : a different GCN adapted for regression with genotype data from the user praxidike97
    - This second GCN implements both classification using the PyTorch geometric [Planetoid](https://pytorch-geometric.readthedocs.io/en/latest/modules/datasets.html#torch_geometric.datasets.Planetoid) dataset and regression using a custom PyTorch geometric dataset of the genotype data.
    - To run on either dataset, uncomment the lines in main, train, and the Net class

Commands to run scripts:
------------------------
- python /multi-omics/GCN/main.py<n>
- python /multi-omics/GCN/praxidike997_gcn_regression.py

If you receive a CUDA blocking error:
-------------------------------------
- CUDA_LAUNCH_BLOCKING=1 python multi-omics/GCN/praxidike97_gcn.py

__[HPCC Development Nodes with GPU cards](https://wiki.hpcc.msu.edu/display/ITH/Development+nodes):__
- dev-intel14-k20	: 20 cores	128GB	Two Nvidia K20 GPUs
- dev-intel16-k80 :	28 cores	256GB	Intel16 node with 4 Nvidia Tesla K80 GPUs
- dev-amd20-v100 : 20 cores	128GB	Two Nvidia K20 GPUs

__Directory of scripts:__ /mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/multi-omics/GCN


__Directory of data:__ /mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic
>>>>>>> 314110fd1adc175794a25433dd6fdad0398f8783
