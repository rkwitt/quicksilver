function dJ = jacobianDetVField(v)
hsize = size(v);
eyeH = eyeHField(hsize(2:end),class(v));
dJ = jacobianDetHField(v+eyeH);

