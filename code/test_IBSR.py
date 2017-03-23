#example file used to predict momentum for IBSR dataset
from collections import Counter
import torch
import torch.nn as nn
from torch.utils.serialization import load_lua
import prediction_network
import util
import numpy as np
import argparse
import h5py


print('config for testing')
print('Loading configuration and network')
config = torch.load('checkpoint10.pth.tar');
patch_size = config['patch_size']
network_feature = config['network_feature']

batch_size = 48

use_multiGPU = True;

n_GPUs = 8;

print('creating net')
net_single = prediction_network.net(network_feature).cuda();
net_single.load_state_dict(config['state_dict'])

if use_multiGPU:
    print('multi GPU net')
    device_ids=range(0, n_GPUs)
    net = torch.nn.DataParallel(net_single, device_ids=device_ids).cuda()
else:
    net = net_single

net.train()

input_batch = torch.zeros(batch_size, 2, patch_size, patch_size, patch_size).cuda()

base_idx = 1;
time_all = 0;

moving_image_dataset = load_lua('IBSR.t7')

for from_image in range(0, moving_image_dataset.size()[0]) :
    for to_image in range(0, target_image_dataset.size()[0]) :
        image_from = moving_image_dataset[from_image].squeeze()
        image_to = moving_image_dataset[to_image].squeeze()
        predict_result = util.predict_momentum(image_from, image_to, input_batch, batch_size, patch_size, net);
        predict_result = predict_result.numpy();
        f = h5py.File("/IBSR18/"+str(from_image+1)+"_"+str(to_image+1)+".h5", "w") 
        dset = f.create_dataset("dataset", data=predict_result)
        f.close()
    #endfor
#endfor






