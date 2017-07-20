# quicksilver
This repository contains code and data for the following paper [[arXiv]](https://arxiv.org/abs/1703.10908):

```
@misc{YangKwittStynerNiethammer17a,
    title        = {Quicksilver: Fast Predictive Image Registration - a Deep Learning Approach},
    author       = {X. Yang and R. Kwitt and M. Styner and M. Niethammer},
    year         = {2017},
    howpublished = {arXiv:1703.10908}}
```


## Disclaimer
**This software is published for academic and non-commercial use only (Apache-2.0)**

## Setup
This code is based or PyTorch and PyCA. It has been tested on Ubuntu 14.04/16.04 LTS with Python 2.7 (using Nvidia TitanX GPUs, with CUDA 8.0).

**Dependencies**:
* [PyTorch](http://pytorch.org/)
* [PyCA](https://bitbucket.org/scicompanat/pyca).
* [CUDA](https://developer.nvidia.com/cuda-downloads)
* [NiftyReg](https://sourceforge.net/projects/niftyreg/)
* [scikit-image](http://scikit-image.org)

*Remark*: NiftyReg is used for optional affine pre-alignment to the ICBM152 atlas before performing deformable registration. scikit-image is used for histogram equalization of the input images. Note: if you use Anaconda as the development platform, you can install scikit-image as suggested in [the Anaconda page](https://anaconda.org/anaconda/scikit-image).

### Exemplary system configuration

We recommend using [Anaconda](https://www.continuum.io/) as your development platform, as this cleanly separates the Python/PyTorch installation (and all required libraries) from your system libraries. We exemplify such an installation below, assuming that the full installation is done under `/scratch` and CUDA (8.0) is installed under `/usr/local/cuda`.

**Install Anaconda 4.3.1 (Python 2.7)**
```bash
mkdir /scratch
cd /scratch
wget https://repo.continuum.io/archive/Anaconda2-4.3.1-Linux-x86_64.sh
bash Anaconda2-4.3.1-Linux-x86_64.sh
```
Follow the installation instructions. We assume that you installed
Anaconda under `/scratch/anaconda2`.

**Install PyCA**

For our code we use a *specific* version of PyCA, so please use the following command to download the PyCA code:
```bash
cd /scratch
git clone git@bitbucket.org:scicompanat/pyca.git
cd ./pyca
git checkout 9954dd5319efaa0ac5f58977e57acf004ad73ed7
mkdir Build && cd Build
ccmake ...
```
Configure (according to your system settings) and compile PyCA.

**Final configuration**

Finally, we clone the Quicksilver repository as
```bash
cd /scratch
git clone https://github.com/rkwitt/quicksilver.git
```
and create a config file `conda.cfg` with the following content:
```bash
export PATH=$PATH:/usr/local/cuda/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64
export PYTHONPATH=$PYTHONPATH:/scratch/pyca/Build/python_module
export CUDA_HOME=/usr/local/cuda
export PATH="/scratch/anaconda2/bin:$PATH"
```
Source this file with `source conda.cfg` to setup the environment.

## Usage

Below is a simple *quickstart* guide on how to use Quicksilver for
*image-to-image* registration (more to come).

### Test data

Four (pre-aligned) example images (from CUMC) can be downloaded [here](https://drive.google.com/open?id=0BxHF82gaPzgSN0IwMnpXTHNibWc).

### Quick start
```bash
cd /scratch/quicksilver/code/applications
python qs_predict.py --moving-image moving_1.nii moving_2.nii moving_3.nii ...
                     --target-image target_1.nii target_2.nii target_3.nii ...
		             --output-prefix prefix_1 prefix_2 prefix_3 ...
```
Add the ``--use-correct`` option if want to use correction network.

To use multiple GPUs for prediction, add ``--n-GPU N`` option, where ``N`` indicates the number of GPUs to use. For maximum performance, combine this option with ``--batch-size K`` to specify the number of 3D patches as network input, and set ``K`` to be a multiple of ``N``.

### Example

In this example, we want to register `m1.nii` (moving) and `m2.nii` (target) from the provided example data. As the images are already pre-aligned, no
additional affine alignment needs to be done.

```bash
cd /scratch2/quicksilver
tar xvfz CUMC_examples.tar.gz
cd code/applications
python qs_predict.py \
    --moving-image ../CUMC_examples/m1.nii \
    --target-image ../CUMC_examples/m2.nii \
    --output-prefix /tmp/
```
This will generate two files: `/tmp/I1.mhd`, `/tmp/I1.raw` (i.e., `m1.nii` aligned to `m2.nii` in the coordinate system of `m2.nii`).

### Use the probablistic version of the network
Use `quicksilver/code/applications/qs_predict_probablistic.py`. There are two main differences compared to qs_predict.py. First, there is an additional option `--samples` that lets you decide the number of times to sample the network (default is 50). The variance of the deformation is samed as 'phiinv_var.mhd'. Second, there are no correction network for the probabilistic network.

### Training a new network from start to finish
1. Do affine alignment and histogram equalization for the training images using `quicksilver/code/tools/preprocessing/affine_and_histogram_eq.py`. If using the default atlas (the MNI152 atlas as `quicksilver/data/atlas/icbm152.nii`), Make sure the training images are in the same coordinate system.
2. Perform registrations on the training data using ```quicksilver/code/tools/LDDMM_optimization/CAvmMatching.py```
3. Gather the moving images, target images and initial momentum from LDDMM optimization for training the network. This is done by using `quicksilver/code/tools/create_pth.py` to gather the images and momentums into .pth.tar files. For example:
```
cd quicksilver/code/tools/
python ./create_pth.py --files moving_image_1 moving_image_2 moving_image_3 --output moving_image_all.pth.tar
python ./create_pth.py --files target_image_1 target_image_2 target_image_3 --output target_image_all.pth.tar
python ./create_pth.py --files momentum_1 momentum_2 momentum_3 --output momentum_all.pth.tar --momentum
```
You can seperated the images/momentums into several .pth.tar files to make each .pth.tar file size reasonable. For example:
```
# make training data part 1
python ./create_pth.py --files moving_image_1 moving_image_2 moving_image_3 --output moving_image_dataset_1.pth.tar
python ./create_pth.py --files target_image_1 target_image_2 target_image_3 --output target_image_dataset_1.pth.tar
python ./create_pth.py --files momentum_1 momentum_2 momentum_3 --output momentum_dataset_1.pth.tar --momentum

# make training data part 2
python ./create_pth.py --files moving_image_4 moving_image_5 moving_image_6 --output moving_image_dataset_2.pth.tar
python ./create_pth.py --files target_image_4 target_image_5 target_image_6 --output target_image_dataset_2.pth.tar
python ./create_pth.py --files momentum_4 momentum_5 momentum_6 --output momentum_dataset_2.pth.tar --momentum
```

Make sure to make the moving images/target images/momentums have the same order in the .pth.tar files.
4. Train the prediction network using `quicksilver/code/tools/qs_train.py`. An example will be
```
cd quicksilver/code/tools
python qs_train.py \
    --moving-image-dataset moving_image_dataset_1.pth.tar moving_image_dataset_2.pth.tar \
    --target-image-dataset target_image_dataset_1.pth.tar target_image_dataset_2.pth.tar \
    --deformation-parameter momentum_dataset_1.pth.tar momentum_dataset_2.pth.tar \
    --deformation-setting-file ./LDDMM_spec.yaml	\
    --output-name ./prediction_network_parameter.pth.tar

```
Here `LDDMM_spec.yaml` defines the setting for the LDDMM optimization algorithm. This information is stored in the network parameter file, and is used when using quicksilver to formulate LDDMM shooting.
5. Create warp-back target image files and momentum difference (between LDDMM optimization and prediction network) files. This is for training the correction network. This operation is done using `quicksilver/code/tools/prepare_correction_training_data.py`. An example would be (suppose we have the training files as `moving_image_1.pth.tar`/`moving_image_2.pth.tar`, `target_image_dataset_1.pth.tar`/`target_image_dataset_2.pth.tar`, and `momentum_dataset_1.pth.tar`/`momentum_dataset_2.pth.tar`):
```
cd quicksilver/code/tools
python prepare_correction_training_data.py \
    --moving-image-dataset moving_image_dataset_1.pth.tar moving_image_dataset_2.pth.tar \
    --target-image-dataset target_image_dataset_1.pth.tar target_image_dataset_2.pth.tar \
    --deformation-parameter momentum_dataset_1.pth.tar momentum_dataset_2.pth.tar \
	--network-parameter ./prediction_network_parameter.pth.tar \
	--warped-back-target-output warped_target_dataset_1.pth.tar warped_target_dataset_1.pth.tar \
	--momentum-residual momentum_diff_1.pth.tar momentum_diff_2.pth.tar
```
6. Train the correction network. The procedure is the same as step 4, except (as in this example here) change `target_image_dataset_x.pth.tar` to `warped_target_dataset_x.pth.tar` and `momentum_dataset_1.pth.tar` to `momentum_diff_1.pth.tar`, and of course change the output file name.


### evaluate your result on the 4 test datasets
1. Perform prediction on the four test datasets (CUMC12, LPBA40, MGH10, IBSR18) use qs_predict.py
2. Calculate the label overlapping score for each test case using `calculate_CUMC_overlap.m`, `calculate_LPBA_overlap.m`, `calculate_IBSR_overlap.m`, `calculate_MGH_overlap.m` in the `quicksilver/code/tools/evaluate_result/` directory
3. Plot the results and compared to the results in the Neuroimage paper by using `quicksilver/code/tools/evaluate_result/generate_label_overlapping_plot.m`