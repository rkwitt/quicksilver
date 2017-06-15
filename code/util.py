import torch
import numpy as np
import h5py
from torch.autograd import Variable

def calculateIdx1D(length, patch_length, step):
	one_dim_pos = torch.arange(0, length-patch_length+1, step)
	if (length-patch_length) % step != 0:
		one_dim_pos = torch.cat((one_dim_pos, torch.ones(1) * (length-patch_length)))
	return one_dim_pos;


def idx2pos(idx, image_size):
	"""
	Given a flattened idx, return the position in the 3D image space.

	Args:

		idx (int): 					Index into flattened 3D volume
		image_size(list of 3 int): 	Size of 3D volume

	"""
	assert(len(image_size)==3)

	pos_x = idx / (image_size[1] * image_size[2]);
	idx_yz = idx % (image_size[1] * image_size[2]);
	pos_y = idx_yz / image_size[2];
	pos_z = idx_yz % image_size[2];
	return torch.LongTensor([pos_x, pos_y, pos_z]);


def pos2idx(pos, image_size):
	"""
	Given a position in the 3D image space, return a flattened idx.

	Args:

		pos (list of 3 int):		Position in 3D volume
		image_size (list of 3 int):	Size of 3D volume
	"""
	assert(len(pos)==3)
	assert(len(image_size)==3)

	return (pos[0] * image_size[1]*image_size[2]) + (pos[1] * image_size[2]) + pos[2];


#given a flatterned idx for a 4D data (n_images * 3D image), return the position in the 4D space
def idx2pos_4D(idx, image_size):
	image_slice = idx / (image_size[0] * image_size[1] * image_size[2])
	single_image_idx = idx % (image_size[0] * image_size[1] * image_size[2])
	single_image_pos = idx2pos(single_image_idx, image_size)
	return torch.cat((image_slice * torch.ones(1).long(), single_image_pos))


# calculate the idx of the patches for 3D dataset (n_images * 3D image)
def calculatePatchIdx3D(num_image, patch_size, image_size, step_size):
	#calculate the idx for 1 3D image
	pos_idx = [calculateIdx1D(image_size[i], patch_size[i], step_size[i]).long() for i in range(0, 3)];
	pos_idx_flat = torch.zeros(pos_idx[0].size()[0] * pos_idx[1].size()[0] * pos_idx[2].size()[0]).long()
	flat_idx = 0;
	pos_3d = torch.zeros(3).long();
	for x_pos in range(0, pos_idx[0].size()[0]):
		for y_pos in range(0, pos_idx[1].size()[0]):
			for z_pos in range(0, pos_idx[2].size()[0]):
				pos_3d[0] = pos_idx[0][x_pos]
				pos_3d[1] = pos_idx[1][y_pos]
				pos_3d[2] = pos_idx[2][z_pos]
				pos_idx_flat[flat_idx] = pos2idx(pos_3d, image_size)
				flat_idx = flat_idx+1;
	
	pos_idx_flat_all = pos_idx_flat.long();

	# calculate the idx across all 3D images in the dataset
	for i in range(1, num_image):
		pos_idx_flat_all = torch.cat((pos_idx_flat_all, pos_idx_flat.long() + i * (image_size[0] * image_size[1] * image_size[2])));

	return pos_idx_flat_all;


# read HDF5 format file
def readHDF5(filename):
	f = h5py.File(filename, 'r')
	data = f['/dataset'][()]
	data = torch.from_numpy(data)
	f.close()
	return data


#convert the image or momentum to prediction space
def convert_to_predict_space(image):
	if (len(image.shape) == 3): #image
		output = np.transpose(image, [2, 1, 0])
	elif (len(image.shape) == 4): # momentum
		output = np.transpose(image, [3, 2, 1, 0])
	else :
		print('does not support 2D yet!')
		sys.exit(1)
	return output;
#enddef


def convert_to_registration_space(image):
	return convert_to_predict_space(image);
#enddef


# predict the momentum given a moving and target image
def predict_momentum(moving, target, input_batch, batch_size, patch_size, net, step_size=14):
    moving = convert_to_predict_space(moving);
    target = convert_to_predict_space(target);

    moving = torch.from_numpy(moving).cuda();
    target = torch.from_numpy(target).cuda();
    data_size = moving.size();
    flat_idx = calculatePatchIdx3D(1, patch_size*torch.ones(3), data_size, step_size*torch.ones(3));
    flat_idx_select = torch.zeros(flat_idx.size());
    #remove the background patches
    for patch_idx in range(1, flat_idx.size()[0]):
        patch_pos = idx2pos_4D(flat_idx[patch_idx], data_size)
        moving_patch = moving[patch_pos[1]:patch_pos[1]+patch_size, patch_pos[2]:patch_pos[2]+patch_size, patch_pos[3]:patch_pos[3]+patch_size]
        target_patch = target[patch_pos[1]:patch_pos[1]+patch_size, patch_pos[2]:patch_pos[2]+patch_size, patch_pos[3]:patch_pos[3]+patch_size]
        if (torch.sum(moving_patch) + torch.sum(target_patch) != 0):
            flat_idx_select[patch_idx] = 1;

    flat_idx_select = flat_idx_select.byte();
    flat_idx = torch.masked_select(flat_idx, flat_idx_select);	

    momentum_predict = torch.zeros(3, data_size[0], data_size[1], data_size[2]).cuda()
    momentum_weight = torch.zeros(3, data_size[0], data_size[1], data_size[2]).cuda()
    #start prediction
    batch_idx = 0;
    while(batch_idx < flat_idx.size()[0]):
        if (batch_idx + batch_size < flat_idx.size()[0]):
            cur_batch_size = batch_size;
        else:
            cur_batch_size = flat_idx.size()[0] - batch_idx

        for slices in range(0, cur_batch_size):
            patch_pos = idx2pos_4D(flat_idx[batch_idx+slices], data_size)
            input_batch[slices, 0] = moving[patch_pos[1]:patch_pos[1]+patch_size, patch_pos[2]:patch_pos[2]+patch_size, patch_pos[3]:patch_pos[3]+patch_size]
            input_batch[slices, 1] = target[patch_pos[1]:patch_pos[1]+patch_size, patch_pos[2]:patch_pos[2]+patch_size, patch_pos[3]:patch_pos[3]+patch_size]

        input_batch_variable = Variable(input_batch, volatile=True)
        recon_batch_variable = net(input_batch_variable)
        for slices in range(0, cur_batch_size):
            patch_pos = idx2pos_4D(flat_idx[batch_idx+slices], data_size)
            momentum_predict[:, patch_pos[1]:patch_pos[1]+patch_size, patch_pos[2]:patch_pos[2]+patch_size, patch_pos[3]:patch_pos[3]+patch_size] += recon_batch_variable.data[slices]
            momentum_weight[:, patch_pos[1]:patch_pos[1]+patch_size, patch_pos[2]:patch_pos[2]+patch_size, patch_pos[3]:patch_pos[3]+patch_size] += 1

        batch_idx += cur_batch_size

    #remove 0 weight areas
    momentum_weight += (momentum_weight == 0).float()
    momentum_predict = momentum_predict.div(momentum_weight).cpu().numpy()

    return convert_to_registration_space(momentum_predict)
#enddef


