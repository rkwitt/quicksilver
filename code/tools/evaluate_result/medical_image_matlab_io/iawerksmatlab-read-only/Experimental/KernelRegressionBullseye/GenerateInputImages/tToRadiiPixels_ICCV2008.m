function r = tToRadiiPixels_ICCV2008(t, imageDims)

minDim        = min(imageDims(:));
rMultRange    = [0.25 0.75];
rInnerMin     = 1/8;
rOuterMax     = 3/9;
rMiddleMin    = rInnerMin+(rOuterMax-rInnerMin)/3;
rOuterMin     = rInnerMin+(rOuterMax-rInnerMin)*2/3;
innerMinR     = rInnerMin*minDim;
innerMaxR     = rMiddleMin*minDim;
middleMinR    = rMiddleMin*minDim;
middleMaxR    = rOuterMin*minDim;
outerMinR     = rOuterMin*minDim;
outerMaxR     = rOuterMax*minDim;

gs1 = rMultRange(1) + f1(t)*(rMultRange(2)-rMultRange(1));
r(:,1) = innerMinR + gs1*(innerMaxR-innerMinR);

gs2 = rMultRange(1) + f2(t)*(rMultRange(2)-rMultRange(1));
r(:,2) = middleMinR + gs2*(middleMaxR-middleMinR);

gs3 = rMultRange(1) + f3(t)*(rMultRange(2)-rMultRange(1));
r(:,3) = outerMinR + gs3*(outerMaxR-outerMinR);
