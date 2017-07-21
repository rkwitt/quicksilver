# Maintainer: Jacob Hinkle <jacob@sci.utah.edu>
pkgname=pyca
pkgver=0.01
pkgrel=1
pkgdesc="Python for Computational Anatomy"
arch=('x86_64')
url="http://bitbucket.org/scicompanat/pyca"
license=('BSD')
groups=()
depends=()
makedepends=('git' 'python2' 'swig' 'fftw')
optdepends=('insight-toolkit: For ITK file support (needs rebuild)'
            'cuda: For GPU computation ability (needs rebuild)')
provides=(pyca)
conflicts=(pyca-git)
replaces=()
backup=()
options=()
install=
source=("https://bitbucket.org/scicompanat/pyca/get/v$pkgver.tar.bz2")
noextract=()
md5sums=('33ebf56758765d582164e7a787bc4920') #generate with 'makepkg -g'

_tardir='scicompanat-pyca-9507157ccc11'

build() {
  #
  # unpack
  #
  cd $srcdir
  tar -jxf v$pkgver.tar.bz2
  mv $_tardir $srcdir/$pkgname-$pkgver

  #
  # set up the build dir
  #
  rm -rf "$srcdir/$pkgname-$pkgver-build"
  mkdir "$srcdir/$pkgname-$pkgver-build"
  cd "$srcdir/$pkgname-$pkgver-build"

  #
  # BUILD HERE
  #
  cmake -D CMAKE_BUILD_TYPE=Release \
        -D PYTHON_EXECUTABLE=/usr/bin/python2 \
        -D PYTHON_INCLUDE_DIR=/usr/include/python2.7 \
        -D BUILD_SHARED_LIBS=ON \
        -D PYTHON_LIBRARY=/usr/lib/libpython2.7.so \
        -D PYTHON_INSTALL_DIR=/usr/lib/python2.7/site-packages \
        -D CUDA_NVCC_FLAGS_RELEASE=--pre-include\ $srcdir/$pkgname-$pkgver/preinc.h $srcdir/$pkgname-$pkgver
  make
}

package() {
  cd "$srcdir/$pkgname-$pkgver-build"
  make DESTDIR=$pkgdir install
}

# vim:set ts=2 sw=2 et:
