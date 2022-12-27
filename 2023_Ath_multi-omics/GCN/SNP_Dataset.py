"""
Custom PyTorch "In Memory" dataset of bi-allelic SNPs for 
383 Arabidopsis thaliana accessions. 

The SNPs are encoded in binary -1 (homozygous alternate allele)
or 1 (homozygous reference allele).

The traits in this dataset are FT (flowering time), RL (leaf number),
RBN (primary branch number), Length (stem length), Diameter 
(flower diameter), and CL (cauline auxiliary branch number)
"""
import os
import torch
from torch_geometric.data import InMemoryDataset, download_url
import pandas as pd
import datatable as dt
os.chdir("/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/multi-omics/GCN")
import calc_adj

class At383(InMemoryDataset):
    """
    Args:
        root (string): Root directory where dataset is saved
        name (string): Name of the dataset
        transform (bool):
        pre_transform (bool):
        
    """
    def __init__(self, root, transform=None, pre_transform=None):
        
        assert split in ["train", "test"]
        
        super().__init__(root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_dir(self) -> str: # root directory
        return os.path.join(self.root, 'raw')

    @property
    def raw_file_names(self): # list of raw data files in the root directory
        return ["SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv",
                "Phenotype_value_383_common_accessions_2017_Grimm.csv", "test_20perc.txt"]

    @property
    def processed_file_names(self):
        return ['data.pt']

    def process(self):
        
        # Read data into huge `Data` list.
        data_list = [...]

        if self.pre_filter is not None:
            data_list = [data for data in data_list if self.pre_filter(data)]

        if self.pre_transform is not None:
            data_list = [self.pre_transform(data) for data in data_list]

        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])

