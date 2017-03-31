# quicksilver
Code and data for paper "Quicksilver: Fast Predictive Image Registration - a Deep Learning Approach"

## Disclaimer 
**This software is published for academic and non-commercial use only.**

## Setup
This code is based or PyTorch and PyCA. It has been tested on Ubuntu 14.04/16.04 LTS with Python 2.7 (using CUDA 8.0 with one 
Nvidia TitanX).

Dependencies:
* [PyTorch](http://pytorch.org/)  
* [PyCA](https://bitbucket.org/scicompanat/pyca). 

For our code we use a *specific* version of PyCA, so please use the following command to download the code:
```
git clone git@bitbucket.org:scicompanat/pyca.git
cd ./pyca
git checkout 9954dd5319efaa0ac5f58977e57acf004ad73ed7
```

CUDA backend:
* [CUDA](https://developer.nvidia.com/cuda-downloads)

Optional:
* [NiftyReg](https://sourceforge.net/projects/niftyreg/). This is used for optional affine alignment to the ICBM152 atlas before performing deformable registration.

### Recommended setup

We recommend using Anaconda 4.3.0, as this cleanly separates the PyTorch installation (and all required libraries) from your
system libraries. We exemplify suchn an installation below:

```
mkdir /scratch
cd /scratch
wget https://repo.continuum.io/archive/Anaconda2-4.3.1-Linux-x86_64.sh
bash Anaconda2-4.3.1-Linux-x86_64.sh
```

## Usage
### Quick start
```
python qs_predict.py --moving-image moving_1.nii moving_2.nii moving_3.nii ...
                     --target-image target_1.nii target_2.nii target_3.nii ...
		     --output-prefix prefix_1 prefix_2 prefix_3 ...
```
Add ``--use-correct`` option if want to use correction network.
                     
To use multiple GPUs for prediction, add ``--n-GPU N`` option, where ``N`` indicates the number of GPUs to use. For maximum performance, combine this option with ``--batch-size K`` to specify the number of 3D patches as network input, and set ``K`` to be a multiple of ``N``.
