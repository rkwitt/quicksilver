function mu = vFieldPointwiseMean(vFields)
% compute componentwise mean of vectors at each point
mu = mean(vFields,ndims(vFields));

