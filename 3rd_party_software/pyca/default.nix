{ stdenv, pkgs,
cmake, swig2, git, which, doxygen,
boost, fftw, fftwFloat,
python, python27Packages,
itk, cudatoolkit, linuxPackages,
useITK ? true, useCUDA ? true
}:

let
    os = stdenv.lib.optionalString;
    optOnOff = pred: ''${os pred "ON"}${os (!pred) "OFF"}'';
in stdenv.mkDerivation rec {
    name = "PyCA-${os useITK "itk-"}${os useCUDA "cuda-"}${version}";
    version = "0.1";
  
    src = ./.;
  
    enableParallelBuilding = true;
  
    nativeBuildInputs = [ cmake swig2 git which doxygen sitePackages ];
    buildInputs = with pkgs; [ boost fftw fftwFloat
    python python27Packages.numpy]
    ++ stdenv.lib.optional useITK [ itk ]
    ++ stdenv.lib.optional useCUDA [ cudatoolkit linuxPackages.nvidia_x11 ]
    ;
  
    cmakeFlags = [
      ''-DBUILD_SHARED_LIBS=OFF''
      ''-DUSE_ITK=${optOnOff useITK}''
      ''-DUSE_CUDA=${optOnOff useCUDA}''
    ];
    # hokey way to get sitePackages directory inside the build environment
    sitePackages = "${python.libPrefix}/site-packages";
    preConfigure = ''
      #export VERBOSE=1
      cmakeFlags="$cmakeFlags -DPYTHON_INSTALL_DIR=$out/lib/$sitePackages
      -DPYTHON_INCLUDE_PATH=${python}/include/python2.7"
    '';
  
    meta = {
      description = "Python for Computational Anatomy";
      homepage = http://bitbucket.org/scicompanat/pyca;
      license = stdenv.lib.licenses.bsd3;
      platforms = stdenv.lib.platforms.all;
    };
}

