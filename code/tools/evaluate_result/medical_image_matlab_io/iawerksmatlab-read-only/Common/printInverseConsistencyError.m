function printInverseConsistencyError(v,vinv)

% first compute error for v \circ vinv
fprintf('Error for v(vinv(.))\n');
fprintf('Min %g, Max %g, Mean %g, Variance %g, Median %g\n',...
  inverseConsistencyError(v,vinv));

% next compute error for vinv \circ v
fprintf('Error for vinv(v(.))\n');
fprintf('Min %g, Max %g, Mean %g, Variance %g, Median %g\n',...
  inverseConsistencyError(vinv,v));

