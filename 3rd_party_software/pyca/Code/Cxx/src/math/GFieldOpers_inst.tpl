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

def RescaleVector = [
true,
false]

template void GFieldOpers::
ComposeHH<${HOperBG}>(Field3D&, const Field3D&, const Field3D&, StreamT);

template void GFieldOpers::
ComposeVH<${VOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void GFieldOpers::
ComposeVInvH<${VOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void GFieldOpers::
ComposeHV<${HOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void GFieldOpers::
ComposeHVInv<${HOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void GFieldOpers::
ComposeTranslation<${FOperBG}>(Field3D&, const Field3D&, const Vec3Df&, StreamT,bool);

template void GFieldOpers::
ApplyH<${FOperBG}>(Field3D&, const Field3D&, const Field3D&, StreamT st);

template void GFieldOpers::
ApplyV<${FOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool s);

template void GFieldOpers::
ApplyVInv<${FOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool s);

template void GFieldOpers::
Resample<${FOperBG}, ${RescaleVector}>(Field3D&, const Field3D&, StreamT);

template void GFieldOpers::
FixedPointInverse<${HOperBG}>(Field3D &ginv, const Field3D& g, unsigned int numIter, StreamT stream, bool onDev);

template void GFieldOpers::
Ad<${VOperBG}>(Field3D& Z, const Field3D& g, const Field3D& X, StreamT s,bool onDev);

template void GFieldOpers::
CoAd<${VOperBG}>(Field3D& n, const Field3D& g, const Field3D& m, StreamT s,bool onDev);

