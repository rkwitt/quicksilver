# quicksilver
This repository contains code and data for the following paper:

```
@misc{YangKwittNiethammer17a,
    title        = {Quicksilver: Fast Predictive Image Registration - a Deep Learning Approach},
    author       = {X. Yang and R. Kwitt and M. Niethammer},
    year         = {2017},
    howpublished = {arXiv:1703.10908}}
```

## Disclaimer
**This software is published for academic and non-commercial use only.**

## Setup
This code is based or PyTorch and PyCA. It has been tested on Ubuntu 14.04/16.04 LTS with Python 2.7 (using CUDA 8.0 with one
Nvidia TitanX).

Dependencies:
* [PyTorch](http://pytorch.org/)
* [PyCA](https://bitbucket.org/scicompanat/pyca).
* [CUDA](https://developer.nvidia.com/cuda-downloads)
* [NiftyReg](https://sourceforge.net/projects/niftyreg/)

*Remark*: NiftyReg is used for optional affine alignment to the ICBM152 atlas before performing deformable registration.

### Exemplary system configuration

We recommend using [Anaconda](https://www.continuum.io/) as your development platform, as this cleanly separates the Python/PyTorch installation (and all required libraries) from your system libraries. We exemplify such an installation below, assuming that the full installation is done under `/scratch` and CUDA (8.0) is installed under `/usr/local/cuda`.

**Install Anaconda 4.3.1 (Python 2.7)**
```
mkdir /scratch
cd /scratch
wget https://repo.continuum.io/archive/Anaconda2-4.3.1-Linux-x86_64.sh
bash Anaconda2-4.3.1-Linux-x86_64.sh
```
Follow the installation instructions. We assume that you installed
Anaconda under `/scratch/anaconda2`.

**Install PyCA**

For our code we use a *specific* version of PyCA, so please use the following command to download the PyCA code:
```
cd /scratch
git clone git@bitbucket.org:scicompanat/pyca.git
cd ./pyca
git checkout 9954dd5319efaa0ac5f58977e57acf004ad73ed7
mkdir Build && cd Build
ccmake ...
```
Configure and compile PyCA.

**Final configuration**

Finally, we clone the Quicksilver repository as
```
cd /scratch
git clone https://github.com/rkwitt/quicksilver.git
```
and create a config file `conda.cfg` with the following content:
```
export PATH=$PATH:/usr/local/cuda/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64
export PYTHONPATH=$PYTHONPATH:/scratch/pyca/Build/python_module
export CUDA_HOME=/usr/local/cuda
export PATH="/scratch/anaconda2/bin:$PATH"
```
Source this file with `source conda.cfg` to setup the environment.

## Usage
### Quick start
```
cd /scratch/quicksilver/code
python qs_predict.py --moving-image moving_1.nii moving_2.nii moving_3.nii ...
                     --target-image target_1.nii target_2.nii target_3.nii ...
		     --output-prefix prefix_1 prefix_2 prefix_3 ...
```
Add ``--use-correct`` option if want to use correction network.

To use multiple GPUs for prediction, add ``--n-GPU N`` option, where ``N`` indicates the number of GPUs to use. For maximum performance, combine this option with ``--batch-size K`` to specify the number of 3D patches as network input, and set ``K`` to be a multiple of ``N``.
