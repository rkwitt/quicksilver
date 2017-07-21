#!/bin/env python2

# NOTE: this was converted from GenPythonModule.sh and updated to be a little
# more cross-platform. JDH Nov 2013

import shutil

import distutils.dir_util

import glob

import os
import sys

def touch(fname, times=None):
    with file(fname, 'a'):
        os.utime(fname, times)


def find_and_delete_endswith(d, suffices, l=0):
    for root, dirs, files in os.walk(d):
        for f in files:
            for s in suffices:
                if f.endswith(s):
                    os.remove(os.path.join(root, f))
        for dr in dirs:
            find_and_delete_endswith(os.path.join(root,dr), suffices, l=l+1)

def find_and_delete_dirs(d, dirnames):
    for root, dirs, files in os.walk(d):
        for f in files:
            if f in dirnames and os.path.isdir(f):
                shutil.rmtree(f)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'Usage: {0} <src_dir> <build_dir>'.format(sys.argv[0])
        sys.exit(1)

    SRC_DIR = sys.argv[1]
    BUILD_DIR = sys.argv[2]

    BIN_DIR = os.path.join(BUILD_DIR, 'bin')
    PYTHON_MOD_DIRNAME = 'python_module'
    PYTHON_MOD_DIR = os.path.join(BUILD_DIR, PYTHON_MOD_DIRNAME)
    # static python code directory:
    PYTHON_SRC_DIR = os.path.join(SRC_DIR, 'Code', 'Python')
    # generated python code directory:
    PYTHON_GEN_DIR = os.path.join(BUILD_DIR, 'Python')
    if not os.path.isdir(PYTHON_MOD_DIR):
        os.mkdir(PYTHON_MOD_DIR)

    touch(os.path.join(PYTHON_MOD_DIR, 'timestamp'))

    MODULEDIR = os.path.join(PYTHON_MOD_DIR, 'PyCA')

    if os.path.isdir(MODULEDIR):
        shutil.rmtree(MODULEDIR, ignore_errors=True)

    distutils.dir_util.copy_tree(PYTHON_SRC_DIR, MODULEDIR)
    distutils.dir_util.copy_tree(PYTHON_GEN_DIR, MODULEDIR)

    # remove .svn dirs
    find_and_delete_dirs('.', ['.svn'])
    find_and_delete_endswith('.', ['.in', '~'])

    #
    # Copy in the python wrapper file
    #
    shutil.copyfile(os.path.join(BUILD_DIR, 'SWIG', 'Core.py'), 
            os.path.join(MODULEDIR, 'Core.py'))

    #
    # Copy in the libraries
    #
    os.chdir(BIN_DIR)
    for ext in ['so','pyd','dll','dylib']:
        for f in glob.glob('_*.'+ext):
            shutil.copyfile(f, os.path.join(MODULEDIR, f))

