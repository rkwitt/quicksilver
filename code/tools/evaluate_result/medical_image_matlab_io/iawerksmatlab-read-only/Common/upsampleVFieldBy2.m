function vOut = upsampleVFieldBy2(vIn)

vOut = zeros([size(vIn,1) size(vIn,2)*2 size(vIn,3)*2],'single');
vOut(:,1:2:end,1:2:end) = vIn;
vOut(:,2:2:end,2:2:end) = vIn;
vOut(:,1:2:end,2:2:end) = vIn;
vOut(:,2:2:end,1:2:end) = vIn;
vOut = vOut * 2;