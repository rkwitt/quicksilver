
def HOperBG = [
BACKGROUND_STRATEGY_CLAMP,
BACKGROUND_STRATEGY_ID,
BACKGROUND_STRATEGY_PARTIAL_ID]

def BoolArg = [true, false]

template void Opers::ComposeHTranslation<${HOperBG}>(Field3D&, const Field3D&, const Vec3Df&, StreamT, bool);

template void Opers::
ResampleH<${HOperBG}, ${BoolArg}>(Field3D&, const Field3D&, StreamT);
