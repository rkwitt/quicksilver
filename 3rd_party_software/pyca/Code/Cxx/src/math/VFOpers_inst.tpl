def VOperBG = [
BACKGROUND_STRATEGY_CLAMP,
BACKGROUND_STRATEGY_WRAP,
BACKGROUND_STRATEGY_ZERO,
BACKGROUND_STRATEGY_PARTIAL_ZERO]

def BoolArg = [true, false]

template void Opers::
ComposeVTranslation<${VOperBG}>(Field3D& h, const Field3D& f, const Vec3Df& t, StreamT stream, bool onDev);

template void Opers::
ResampleV<${VOperBG}, ${BoolArg}>(Field3D&, const Field3D&, StreamT);
