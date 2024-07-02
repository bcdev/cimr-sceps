%-------------------------------------------------------------------------------
%
% SUBMODULE   A function to add random noise to
%	      a mean state fullfilling its full
%	      covariance matrix.
%
% FORMAT   x = rnd_generator( mean_x, S_x, no_cases )   
%        
% OUT   x	  Randomly generated states  
%
% IN    mean_x    Mean state		
%	S_x	  Full covariance matrix 
%	no_cases  Number of cases to generate
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-01-02
%-------------------------------------------------------------------------------

function x = rnd_generator( mean_x, S_x, no_cases )


if isdiag(S_x)

  %= simple method

  randn(no_cases,1);
  dS_x = diag( S_x );
  ns   = length(dS_x);
  x    = zeros(no_cases,ns);

  for f = 1:ns

    x(:,f) = mean_x(f) + sqrt(dS_x(f)) * randn(no_cases,1);

  end

else

  %= full cov matrix, choleski decomposition
 
 [R,p] = chol( S_x );

  %= making cov matrix positive definite

  if p ~= 0

    [V,D] = eig( cov_emis );     % Calculate the eigendecomposition 
				 % of the matrix (A = V*D*V') 
  	                         % where "D" is a diagonal matrix 
				 % holding the eigenvalues of your 
				 % matrix "A"

     d= diag(D);       		 % Get the eigenvalues in a vector "d" 

     d(d <= 1e-7) = 1e-7;  	 % Set any eigenvalues that are lower 
				 % than threshold "TH" ("TH" here being 
                        	 % equal to 1e-7) to a fixed non-zero 
				 % "small" value (here assumed equal 
				 % to 1e-7)

     D_c = diag(d);        	 % Built the "corrected" diagonal matrix 
				 % "D_c"

     cov_emis = V*D_c*V';      	 % Recalculate your matrix "A" in its 
				 % PD variant "A_PD"

  end


  %=== rnd generation of no_cases

  n = size(S_x,1);
  X = randn( no_cases, n );
  X = X - mean(X);
  X = X * inv(chol(cov(X)));
  X = X * chol(S_x); 

  x  = repmat(mea_x', r, 1);
  x   = mea + X;

end


return 
