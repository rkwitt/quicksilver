function v = greensFunction(b,method,parameters)

if isMatlab
  stderr = 2;
end

%method = 'fluid';
%parameters = [0.1 0 1];

switch ndims(b)-1
 case 3,
  % compute dft
  fftx = fftn(squeeze(b(1,:,:,:)));
  ffty = fftn(squeeze(b(2,:,:,:)));
  fftz = fftn(squeeze(b(3,:,:,:)));
  L = size(b,2);
  M = size(b,3);
  N = size(b,4);
  switch method
   case 'fluid'
    fprintf('F');
    alpha = parameters(1);
    beta = parameters(2);
    gamma = parameters(3);

    % apply operator to dft coeffs
    for u = 0:L-1
      for v = 0:M-1
        for w = 0:N-1
	  Lambda(1,1) = ...
	     -2*(alpha + beta)*cos(2*pi*u/L)...
	     -2*alpha*cos(2*pi*v/M)...
	     -2*alpha*cos(2*pi*w/N)...
	     +6*alpha+2*beta+gamma;
	  Lambda(2,2) = ...
	     -2*alpha*cos(2*pi*u/L)...
	     -2*(alpha + beta)*cos(2*pi*v/M)...
	     -2*alpha*cos(2*pi*w/N)...
	     +6*alpha+2*beta+gamma;
	  Lambda(3,3) = ...
	     -2*alpha*cos(2*pi*u/L)...
	     -2*alpha*cos(2*pi*v/M)...
	     -2*(alpha + beta)*cos(2*pi*w/N)...
	     +6*alpha+2*beta+gamma;
	  Lambda(1,2) = beta*sin(2*pi*u/L)*sin(2*pi*v/M);
	  Lambda(2,1) = Lambda(1,2);
	  Lambda(1,3) = beta*sin(2*pi*u/L)*sin(2*pi*w/N);
	  Lambda(3,1) = Lambda(1,3);
	  Lambda(2,3) = beta*sin(2*pi*v/M)*sin(2*pi*w/N);
	  Lambda(3,2) = Lambda(2,3);


	  % solve for v = Lambda^{-1} b
	  bx = [fftx(u+1,v+1,w+1); ffty(u+1,v+1,w+1); fftz(u+1,v+1,w+1)];
	  vx = Lambda \ bx;
	  fftx(u+1,v+1,w+1) = vx(1);
	  ffty(u+1,v+1,w+1) = vx(2);
	  fftz(u+1,v+1,w+1) = vx(3);
        end
      end     
    end
   case 'laplace'
    fprintf('L');
    alpha = parameters(1);
    beta  = parameters(2);
    gamma = parameters(3);
    % apply operator to dft coeffs
    for u = 0:L-1
      for v = 0:M-1
        for w = 0:N-1
	  Lambda = eye(3);
	  Lambda(1,1) = ...
	     -2*alpha*(cos(2*pi*u/L)+cos(2*pi*v/M)+cos(2*pi*w/N))...
	     +6*alpha+gamma;
	  Lambda(2,2) = ...
	     -2*alpha*(cos(2*pi*u/L)+cos(2*pi*v/M)+cos(2*pi*w/N))...
	     +6*alpha+gamma;
	  Lambda(3,3) = ...
	     -2*alpha*(cos(2*pi*u/L)+cos(2*pi*v/M)+cos(2*pi*w/N))...
	     +6*alpha+gamma;

	  % solve for v = Lambda^{-1} b
	  bx = [fftx(u+1,v+1,w+1); ffty(u+1,v+1,w+1); fftz(u+1,v+1,w+1)];
	  vx = Lambda \ bx;
	  fftx(u+1,v+1,w+1) = vx(1);
	  ffty(u+1,v+1,w+1) = vx(2);
	  fftz(u+1,v+1,w+1) = vx(3);
        end
      end     
    end
   case 'gaussian'
   end
  % now take inverse dft of reslulting coefficients
  v = zeros(size(b), class(b));
  v(1,:,:,:) = real(ifftn(fftx));
  v(2,:,:,:) = real(ifftn(ffty));
  v(3,:,:,:) = real(ifftn(fftz));
 case 2,
  % compute dft
  fftx = fftn(squeeze(b(1,:,:)));
  ffty = fftn(squeeze(b(2,:,:)));
  L = size(b,2);
  M = size(b,3);
  switch method
   case 'fluid'
    fprintf('F');
    alpha = parameters(1);
    beta = parameters(2);
    gamma = parameters(3);

    % apply operator to dft coeffs
    for u = 0:L-1
      for v = 0:M-1
	Lambda(1,1) = ...
	     -2*(alpha + beta)*cos(2*pi*u/L)...
	     -2*alpha*cos(2*pi*v/M)...
	     +4*alpha+2*beta+gamma;
	Lambda(2,2) = ...
	     -2*alpha*cos(2*pi*u/L)...
	     -2*(alpha + beta)*cos(2*pi*v/M)...
	     +4*alpha+2*beta+gamma;
	Lambda(1,2) = beta*sin(2*pi*u/L)*sin(2*pi*v/M);
	Lambda(2,1) = Lambda(1,2);

	% solve for v = Lambda^{-1} b
	bx = [fftx(u+1,v+1); ffty(u+1,v+1)];
	vx = Lambda \ bx;
	fftx(u+1,v+1) = vx(1);
	ffty(u+1,v+1) = vx(2);
      end     
    end
   case 'laplace'
    fprintf(stderr,'L');
    alpha = parameters(1);
    beta  = parameters(2);
    gamma = parameters(3);
    % apply operator to dft coeffs
    for u = 0:L-1
      for v = 0:M-1
	Lambda = eye(2);
	Lambda(1,1) = ...
	   -2*alpha*(cos(2*pi*u/L)+cos(2*pi*v/M))...
	   +4*alpha+gamma;
	Lambda(2,2) = ...
	   -2*alpha*(cos(2*pi*u/L)+cos(2*pi*v/M))...
	   +4*alpha+gamma;

	% solve for v = Lambda^{-1} b
	bx = [fftx(u+1,v+1); ffty(u+1,v+1)];
	vx = Lambda \ bx;
	fftx(u+1,v+1) = vx(1);
	ffty(u+1,v+1) = vx(2);
      end     
    end
   case 'gaussian'
   end
  % now take inverse dft of reslulting coefficients
  v = zeros(size(b), class(b));
  v(1,:,:) = real(ifftn(fftx));
  v(2,:,:) = real(ifftn(ffty));
 otherwise,
  error('im must be a 2 or 3 dimensional array');
end 


