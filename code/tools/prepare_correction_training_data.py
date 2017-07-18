# add LDDMM shooting code into path
import sys
sys.path.append('../vectormomentum/Code/Python');
sys.path.append('../library')

from subprocess import call
import argparse
import os.path

#Add deep learning related libraries
from collections import Counter
import torch
import prediction_network
import util
import numpy as np
from skimage import exposure

#Add LDDMM registration related libraries
# pyca modules
import PyCA.Core as ca
import PyCA.Common as common
#import PyCA.Display as display
# vector momentum modules
# others
import logging
import copy
import math
import gc
import registration_methods

parser = argparse.ArgumentParser(description='Use trained prediction network to generate data for training the correction network')
requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument('--moving-image-dataset', nargs='+', required=True, metavar=('m1', 'm2, m3...'),
                           help='List of moving images datasets stored in .pth.tar format (or .t7 format for the old experiments in the Neuroimage paper). File names are seperated by space.')
requiredNamed.add_argument('--target-image-dataset', nargs='+', required=True, metavar=('t1', 't2, t3...'),
                           help='List of target images datasets stored in .pth.tar format (or .t7 format for the old experiments in the Neuroimage paper). File names are seperated by space.')
requiredNamed.add_argument('--deformation-parameter', nargs='+', required=True, metavar=('o1', 'o2, o3...'),
                           help='List of target deformation parameter files to predict to, stored in .pth.tar format (or .t7 format for the old experiments in the Neuroimage paper). File names are seperated by space.')
requiredNamed.add_argument('--network-parameter', required=True, metavar=('parameter.pth.tar'),
                           help='network parameter for the prediction network')
requiredNamed.add_argument('--warped-back-target-output', nargs='+', required=True, metavar = ('w1', 'w2, w3...'),
                           help='directory+names of warped back target image datasets, in .pth.tar format. The datapart seperation is the same as the original training image dataparts.')
requiredNamed.add_argument('--momentum-residual', nargs='+', required=True, metavar = ('r1', 'r2, r3...'),
                           help='directory+names of residual of momentum between optimization and prediction, in .pth.tar format. The datapart seperation is the same as the original training momentum dataparts.')

parser.add_argument('--batch-size', type=int, default=64, metavar='N',
                    help='input batch size for prediction network (default: 64)')
parser.add_argument('--n-GPU', type=int, default=1, metavar='N',
                    help='number of GPUs used for prediction (default: 1). For maximum efficiency please set the batch size divisible by the number of GPUs.')
parser.add_argument('--shoot-steps', type=int, default=0, metavar='N',
                    help='time steps for geodesic shooting. Ignore this option to use the default step size used by the registration model.')

args = parser.parse_args()


def check_args(args):
    len_files = len(args.moving_image_dataset)
    if len_files != len(args.target_image_dataset):
        print('The number of the target image dataset files is not consistent with the moving image dataset files.')
    if len_files != len(args.deformation_parameter):
        print('The number of the deformation parameter dataset files is not consistent with the moving image dataset files.')
    if len_files != len(args.warped_back_target_output):
        print('The number of the warp back target image dataset files is not consistent with the moving image dataset files.')
    if len_files != len(args.momentum_residual):
        print('The number of the momentum residual dataset files is not consistent with the moving image dataset files.')


def create_net(args, network_config):
    net_single = prediction_network.net(network_config['network_feature']).cuda();
    net_single.load_state_dict(network_config['state_dict'])

    if (args.n_GPU > 1) :
        device_ids=range(0, args.n_GPU)
        net = torch.nn.DataParallel(net_single, device_ids=device_ids).cuda()
    else:
        net = net_single

    net.train()
    return net;
#enddef

def predict_dataset(args):
    #create prediction network
    network_config = torch.load(args.network_parameter)
    net = create_net(args, network_config)
    
    predict_transform_space = False
    if 'matlab_t7' in network_config:
        predict_transform_space = True

    batch_size = args.batch_size
    patch_size = network_config['patch_size']
    input_batch = torch.zeros(batch_size, 2, patch_size, patch_size, patch_size).cuda()
    for datapart_idx in range(0, len(args.moving_image_dataset)):
        predict_each_datapart(args, net, network_config, input_batch, datapart_idx, batch_size, patch_size, predict_transform_space)
        gc.collect()
#enddef

def predict_each_datapart(args, net, network_config, input_batch, datapart_idx, batch_size, patch_size, predict_transform_space):
    moving_image = torch.load(args.moving_image_dataset[datapart_idx])
    target_image = torch.load(args.target_image_dataset[datapart_idx])
    optimization_momentum = torch.load(args.deformation_parameter[datapart_idx])
    for slice_idx in range(0, moving_image.size()[0]):
        print(slice_idx)
        moving_slice = moving_image[slice_idx].numpy()
        target_slice = target_image[slice_idx].numpy()
        if predict_transform_space:
            moving_slice = util.convert_to_registration_space(moving_slice)
            target_slice = util.convert_to_registration_space(target_slice)

        predicted_momentum = util.predict_momentum(moving_slice, target_slice, input_batch, batch_size, patch_size, net, predict_transform_space);
        m0_reg = common.FieldFromNPArr(predicted_momentum['image_space'], ca.MEM_DEVICE);
            
        moving_image_ca = common.ImFromNPArr(moving_slice, ca.MEM_DEVICE)
        target_image_ca = common.ImFromNPArr(target_slice, ca.MEM_DEVICE)

        registration_result = registration_methods.geodesic_shooting(moving_image_ca, target_image_ca, m0_reg, args.shoot_steps, ca.MEM_DEVICE, network_config)
            
        target_inv = common.AsNPCopy(registration_result['I1_inv'])
        print(target_inv.shape)
        if predict_transform_space:
            target_inv = util.convert_to_predict_space(target_inv)
        print(target_inv.shape)
        target_inv = torch.from_numpy(target_inv)

        target_image[slice_idx] = target_inv

        optimization_momentum[slice_idx] = optimization_momentum[slice_idx] - torch.from_numpy(predicted_momentum['prediction_space'])
        
    torch.save(target_image, args.warped_back_target_output[datapart_idx])
    torch.save(optimization_momentum, args.momentum_residual[datapart_idx])

if __name__ == '__main__':
    check_args(args)
    predict_dataset(args)