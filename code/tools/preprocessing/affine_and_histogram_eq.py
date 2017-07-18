# add LDDMM shooting code into path
import sys
sys.path.append('../../vectormomentum/Code/Python');
from subprocess import call
import argparse
import os.path

#Add LDDMM registration related libraries
# pyca modules
import PyCA.Core as ca
import PyCA.Common as common

import numpy as np
from skimage import exposure

parser = argparse.ArgumentParser(description='Perform preprocessing (affine transformation, intensity normalization (0 to 1) and histogram equalization (optional)) to input images')

requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument('--input-images', nargs='+', required=True, metavar=('i1', 'i2, i3...'),
                           help='List of input images, seperated by space.')
requiredNamed.add_argument('--output-images', nargs='+', required=True, metavar=('o1', 'o2, o3...'),
                           help='List of output image directory and file names, seperated by space. Do not save as .nii format as PyCA will flip the image. The recommended format to save images is .mhd')
parser.add_argument('--atlas', default="../../../data/atlas/icbm152.nii",
                    help="Atlas to use for (affine) pre-registration")

parser.add_argument('--input-labels', nargs='+', metavar=('l_i1', 'l_i2, l_i3...'),
                           help='List of input label maps for the input images, seperated by space.')

parser.add_argument('--output-labels', nargs='+', metavar=('l_o1', 'l_o2, l_o3...'),
                           help='List of output label maps, seperated by space.')
parser.add_argument('--histeq', action='store_true', default=False,
                    help='Perform histogram equalization to the moving and target images.')

args = parser.parse_args()

def check_args(args):
    if (len(args.input_images) != len(args.output_images)):
        print('The number of input images is not consistent with the number of output images!')
        sys.exit(1)
    if ((args.input_labels is None) ^ (args.output_labels is None)):
        print('The input labels and output labels need to be both defined!')
        sys.exit(1)
    if ((args.input_labels is not None) and (len(args.input_labels) != len(args.output_labels))):
        print('The number of input labels is not consistent with the number of output labels!')

def affine_transformation(args):
    for i in range(0, len(args.input_images)):
        call(["reg_aladin",
            "-noSym", "-speeeeed", "-ref", args.atlas ,
            "-flo", args.input_images[i],
            "-res", args.output_images[i],
            "-aff", args.output_images[i]+'_affine_transform.txt'])     
        if (args.input_labels is not None):
            call(["reg_resample",
                "-ref", args.atlas,
                "-flo", args.input_labels[i],
                "-res", args.output_labels[i],
                "-trans", args.output_images[i]+'_affine_transform.txt',
                "-inter", str(0)]) 

def intensity_normalization_histeq(args):
    for i in range(0, len(args.input_images)):
        image = common.LoadITKImage(args.output_images[i], ca.MEM_HOST)
        grid = image.grid()
        image_np = common.AsNPCopy(image)
        nan_mask = np.isnan(image_np)
        image_np[nan_mask] = 0
        image_np /= np.amax(image_np)

        # perform histogram equalization if needed
        if args.histeq:
            image_np[image_np != 0] = exposure.equalize_hist(image_np[image_np != 0])
        image_result = common.ImFromNPArr(image_np, ca.MEM_HOST);
        image_result.setGrid(grid)
        common.SaveITKImage(image_result, args.output_images[i])

if __name__ == '__main__':
    print((args.input_labels is None) and (args.output_labels is None))
    
    affine_transformation(args);
    intensity_normalization_histeq(args)
