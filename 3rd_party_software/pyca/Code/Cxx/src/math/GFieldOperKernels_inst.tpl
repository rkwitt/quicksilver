def HOperBG = [
BACKGROUND_STRATEGY_ID,
BACKGROUND_STRATEGY_PARTIAL_ID,
BACKGROUND_STRATEGY_CLAMP]

def VOperBG = [
BACKGROUND_STRATEGY_ZERO,
BACKGROUND_STRATEGY_PARTIAL_ZERO,
BACKGROUND_STRATEGY_WRAP,
BACKGROUND_STRATEGY_CLAMP]

def FOperBG = [
BACKGROUND_STRATEGY_ZERO,
BACKGROUND_STRATEGY_PARTIAL_ZERO,
BACKGROUND_STRATEGY_ID,
BACKGROUND_STRATEGY_PARTIAL_ID,
BACKGROUND_STRATEGY_WRAP,
BACKGROUND_STRATEGY_CLAMP]

def Fwd = [
true,
false]

def RescaleVector = [
true,
false]

template
void ComposeVH<${Fwd}, ${VOperBG}>
(float* d_hx, 
 float* d_hy, 
 float* d_hz,
 const float* d_vx, 
 const float* d_vy, 
 const float* d_vz,
 const float* d_gx, 
 const float* d_gy, 
 const float* d_gz,
 const float& delta, 
 int w, int h, int l,
 float spX, float spY, float spZ, 
 StreamT stream, bool onDev);

template
void ComposeHV<${Fwd}, ${HOperBG}>
(float* d_hx, 
 float* d_hy, 
 float* d_hz,
 const float* d_gx, 
 const float* d_gy, 
 const float* d_gz,
 const float* d_vx, 
 const float* d_vy, 
 const float* d_vz,
 const float& delta,
 int w, int h, int l,
 float spX, float spY, float spZ, 
 StreamT stream, bool onDev);

template
void
ComposeTranslation<${FOperBG}>
(float *d_ox, 
 float *d_oy, 
 float *d_oz, 
 const float *d_ix, 
 const float *d_iy, 
 const float *d_iz, 
 const Vec3Di& sz,
 const Vec3Df& t, 
 StreamT stream, 
 bool onDev) ;

template
void
ApplyH<${FOperBG}>
(float *d_ox, 
 float *d_oy, 
 float *d_oz, 
 const float *d_ix, 
 const float *d_iy, 
 const float *d_iz, 
 const float *d_hx, 
 const float *d_hy, 
 const float *d_hz, 
 const Vec3Di &sz,
 StreamT stream);

template
void ApplyV<${Fwd}, ${FOperBG}>
(float *d_ox, 
 float *d_oy, 
 float *d_oz, 
 const float *d_ix, 
 const float *d_iy, 
 const float *d_iz, 
 const float *d_ux, 
 const float *d_uy, 
 const float *d_uz, 
 const Vec3Di &sz,
 const Vec3Df &sp,
 const float& delta, 
 StreamT stream, bool onDev);

template
void
Resample<${FOperBG},  ${RescaleVector}>
(float *d_ox, 
 float *d_oy, 
 float *d_oz, 
 const Vec3Di &oSz,
 const float *d_ix, 
 const float *d_iy, 
 const float *d_iz, 
 const Vec3Di &iSz,
 StreamT stream);

template
void
FixedPointInverse<${HOperBG}> 
(float *ginvx, 
 float *ginvy, 
 float *ginvz, 
 const float *gx, 
 const float *gy, 
 const float *gz, 
 const Vec3Di &sz,
 unsigned int numIter, 
 StreamT stream, bool onDev);

template
void
Ad<${VOperBG}> 
(float *Zx, 
 float *Zy, 
 float *Zz, 
 const float *gx, 
 const float *gy, 
 const float *gz, 
 const float *Xx,
 const float *Xy,
 const float *Xz,
 const Vec3Di &sz,
 const Vec3Df &sp,
 StreamT s,bool onDev);

template
void
CoAd<${VOperBG}> 
(float *nx, 
 float *ny, 
 float *nz, 
 const float *gx, 
 const float *gy, 
 const float *gz, 
 const float *mx,
 const float *my,
 const float *mz,
 const Vec3Di &sz,
 const Vec3Df &sp,
 StreamT s,bool onDev);
