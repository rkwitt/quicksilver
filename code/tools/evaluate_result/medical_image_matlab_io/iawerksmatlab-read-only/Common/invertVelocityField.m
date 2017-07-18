function vinv = invertVelocityField(v,method)

if nargin < 2
  method = 'zero';
end

switch method
  case 'zero'
    vinv = -v;
  case 'iter'
    switch ndims(v)-1
      case 3
        vsize = [size(v,2) size(v,3) size(v,4)];
        maxIter = 10;
        convergenceThreshold = 0.1;
        id = eyeHField(vsize);
        ax = zeros(vsize,class(v));
        ay = zeros(vsize,class(v));
        az = zeros(vsize,class(v));
        converged = false;
        
        % debugging
        %vinv(1,:,:) = -ax;
        %vinv(2,:,:) = -ay;
        %printInverseConsistencyError(v,vinv);
         vinv(1,:,:,:) = -ax;
         vinv(2,:,:,:) = -ay;
         vinv(3,:,:,:) = -az;
         ice = inverseConsistencyError(v,vinv);
         maxICE = ice(2);
         fprintf(':0[ice=%0.3f] ',maxICE);
         
        for i=1:maxIter
          fprintf('%d',i);
          % save a copy to check for convergence later
          oldax = ax;
          olday = ay;
          oldaz = az;
          
          % iterate
          ax = interp3(squeeze(v(1,:,:,:)),squeeze(id(2,:,:,:))-ay,squeeze(id(1,:,:,:))-ax,squeeze(id(3,:,:,:))-az);
          ay = interp3(squeeze(v(2,:,:,:)),squeeze(id(2,:,:,:))-ay,squeeze(id(1,:,:,:))-ax,squeeze(id(3,:,:,:))-az);
          az = interp3(squeeze(v(3,:,:,:)),squeeze(id(2,:,:,:))-ay,squeeze(id(1,:,:,:))-ax,squeeze(id(3,:,:,:))-az);
          
          % check for convergence
          dx = oldax-ax;
          dy = olday-ay;
          dz = oldaz-az;
          maxChange = sqrt(max(max(max(dx.*dx + dy.*dy + dz.*dz))));
          
          vinv(1,:,:,:) = -ax;
          vinv(2,:,:,:) = -ay;
          vinv(3,:,:,:) = -az;
          ice = inverseConsistencyError(v,vinv);
          maxICE = ice(2);
          fprintf('[ice=%0.3f] ',maxICE);
          
          if maxChange < convergenceThreshold
            converged = true;
            break
          end
        end
        if converged
          fprintf('Converged after %d iterations\n',i);
        else
          fprintf('Finished %d iterations without converging\n',i);
          warning('Finished iterations without converging.');
        end
        
        vinv(1,:,:,:) = -ax;
        vinv(2,:,:,:) = -ay;
        vinv(3,:,:,:) = -az;
      case 2
        vsize = [size(v,2) size(v,3)];
        maxIter = 10;
        convergenceThreshold = 0.1;
        id = eyeHField(vsize);
        ax = zeros(vsize,class(v));
        ay = zeros(vsize,class(v));
        converged = false;

        % debugging
        %vinv(1,:,:) = -ax;
        %vinv(2,:,:) = -ay;
        %printInverseConsistencyError(v,vinv);
         vinv(1,:,:) = -ax;
         vinv(2,:,:) = -ay;
         ice = inverseConsistencyError(v,vinv);
         maxICE = ice(2);
         fprintf(':0[ice=%0.3f] ',maxICE);
        
        for i=1:maxIter
          fprintf('%d',i);
          % save a copy to check for convergence later
          oldax = ax;
          olday = ay;
          
          % iterate
          ax = interp2(squeeze(v(1,:,:)),squeeze(id(2,:,:))-ay,squeeze(id(1,:,:))-ax);
          ay = interp2(squeeze(v(2,:,:)),squeeze(id(2,:,:))-ay,squeeze(id(1,:,:))-ax);
          
          % check for convergence
          dx = oldax-ax;
          dy = olday-ay;
          maxChange = sqrt(max(max(dx.*dx + dy.*dy)));

          vinv(1,:,:) = -ax;
          vinv(2,:,:) = -ay;
          ice = inverseConsistencyError(v,vinv);
          maxICE = ice(2);
          fprintf('[ice=%0.3f] ',maxICE);
          
          if maxChange < convergenceThreshold
            converged = true;
            break
          end
        end
        if converged
          fprintf('Converged after %d iterations\n',i);
        else
          fprintf('Finished %d iterations without converging\n',i);          
          warning('Finished iterations without converging.');
        end
        
        vinv(1,:,:) = -ax;
        vinv(2,:,:) = -ay;
      otherwise
        error('velocity field must be 2 or 3 dimensions');
    end
  otherwise
    error('method must be ''zero'', ''first'', ''iter''');
end

