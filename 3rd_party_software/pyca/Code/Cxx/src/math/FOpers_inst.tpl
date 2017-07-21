def VOperBG = [
BACKGROUND_STRATEGY_CLAMP,
BACKGROUND_STRATEGY_WRAP,
BACKGROUND_STRATEGY_ZERO,
BACKGROUND_STRATEGY_PARTIAL_ZERO]

def HOperBG = [
BACKGROUND_STRATEGY_CLAMP,
BACKGROUND_STRATEGY_ID,
BACKGROUND_STRATEGY_PARTIAL_ID]

def AllBG = [
BACKGROUND_STRATEGY_CLAMP,
BACKGROUND_STRATEGY_WRAP,
BACKGROUND_STRATEGY_ZERO,
BACKGROUND_STRATEGY_PARTIAL_ZERO,
BACKGROUND_STRATEGY_ID,
BACKGROUND_STRATEGY_PARTIAL_ID]

def IDBG = [
BACKGROUND_STRATEGY_ID,
BACKGROUND_STRATEGY_PARTIAL_ID]

def BoolArg = [true, false]

template void Opers::
ComposeVH<${VOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void Opers::
ComposeVInvH<${VOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void Opers::
ComposeHV<${HOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void Opers::
ComposeTranslation<${AllBG}>(Field3D&, const Field3D&, const Vec3Df&, StreamT,bool);

template void Opers::
ApplyH<${AllBG}>(Field3D&, const Field3D&, const Field3D&, StreamT st);

template void Opers::
ApplyV<${AllBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool s);

template void Opers::
ApplyVInv<${AllBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool s);

template void Opers::
Resample<${AllBG}, ${BoolArg}>(Field3D&, const Field3D&, StreamT);

template void Opers::
FixedPointInverse<${IDBG}>(Field3D &ginv, const Field3D& g, unsigned int numIter, StreamT stream, bool onDev);

template void Opers::
Ad<${VOperBG}>(Field3D& Z, const Field3D& g, const Field3D& X, StreamT st, bool onDev);

template void Opers::
CoAd<${VOperBG}>(Field3D& n, const Field3D& g, const Field3D& m, StreamT st, bool onDev);

