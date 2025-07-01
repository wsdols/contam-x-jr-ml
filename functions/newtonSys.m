function [varargout] = newtonSys(Jfun, sim_params, x0, varargin)
% newtonSys  Newton's method for systems of nonlinear equations.
%
% Synopsis:  
%    [x, Zones, Paths] = newtonSys(Jfun,sim_params, x0, arg1, arg2, ...)
%
% Input:  
%   Jfun = (string) name of m-file that returns matrix J and vector F
%   x0   = initial guess at solution vector, x
%          arg1, arg2 = arguments that are passed through to the m-file 
%          defined by 'Jfun', e.g., Jf_NZones_Powerlaw()
%   arg1 = Paths struct, arg2 = Zones struct
%
% Output:  
%   arg1 = x = solution vector;  
%       x is returned after k iterations if tolerances are met, or after maxit
%       iterations if tolerances are not met.
%   arg2 = Paths struct
%   arg3 = Zones struct
%
%--------------------------------------------------------------------------
%                                                                           
% Solve non-linear airflow equations using Newton-Raphson with simple under 
%   relaxation and full matrix LU decomposition w/ no pivoting to solve the 
%   linear equations for each NR iteration. Calls Jfun(), e.g., 
%   Jf_NZones_Powerlaw() to fill the Jacobian Matrix and Residiual 
%   function Vector.           
%                                                                           
% BACKGROUND                                                                
% ----------                                                                
% Dennis and Schnabel. 1996. Numerical Methods for Unconstrained            
% Optimization and Nonlinear Equations.                                     
%   J(Xc) Sn = -F(Xc) => [A]{x} = {b}                                       
%   {Xc} is the current iterate (or initial guess). In CONTAM, this is the  
%        vector of Zone Reference Pressures being determined via this non-  
%        linear solution method.                                            
%   {F} is the the sum of flows out of each node n - this is the function   
%       often referred to as the residual vector being made zero in NR      
%       method. NOTE: CONTAM sign convention for dPnm is such that positive 
%       Flows occur when Pn > Pm and flows are out of Zone n. Leading to    
%       the sign convention for this residual function (see WSD NOTES).     
%       Fn = SUMi(w_nm), where i includes all paths connecting               
%            zone n and m.                                                  
%   [J(Xc)] is an N x N Jacobian matrix                                     
%   {Sn} is the newton iteration obtained by solving the linear system      
%         [J]{Sn} = -{F}                                                    
%   Solve [J]{Sn} = -{F} for Sn and obtain the next iterate X+ from         
%     X+ = Xc + Sn                                                          
%   Check for convergence:                                                  
%    ||X+|| < EPS                                                           
%    NOTE: ContamX uses other methods to test for convergence.              
%                                                                           
%   WSD NOTES:                                                              
%   * D&S formulation for iterate k to determine k+1:                       
%                           [J]{Sk}=-{F} => Xk+1 = Xk + r Sk                
%     Alternate:           -[J]{Sk}= {F} => Xk+1 = Xk + r Sk                
%   * ContamX formulation:  [J]{Sk}= {F} => Xk+1 = Xk - r Sk                
%   * "r" is the relaxation coefficient: 0.0 < r <= 1.0                     
%     GNW found 0.75 works very well for CONTAM problems.                   
%   * Set SYSTEM_SIGN to demonstrate that either D&S or ContamX is correct  
%     as long as consistency is maintained between newton_sys() and         
%     fill_Jf():                                                            
%   * ContamXJr formulation:                                                
%       [J(Xk)]{Sk} = SYSTEM_SIGN*{F(Xk)}   =>                              
%       Xk + 1 = Xk - (SYSTEM_SIGN * r Sk)                                  
%       SYSTEM_SIGN = -1 for ContamXJr, +1 for ContamX                      
%   JACOBIAN ELEMENTS                                                       
%   * DP = Pn - Pm (assuming Ps and Pw are incorporated into Pn and Pm)     
%   * Mass airflow through path i: wi = w_nm = C(Pn - Pm)^x                 
%     w = C y^x                                                             
%     dw/dPn = xC y^(x-1) dy/DP                                             
%            = xC (y^x)/x                                                   
%            =  C x w/DP                                                    
%       DP+ => dw/dPn = +/+ = +                                             
%       DP- => dw/dPn = -/- = +                                             
%     dw/dPm = -C x w/DP                                                    
%       DP+ => dw/dPm = -(+/+) = -                                          
%       DP- => dw/dPm = -(-/-) = -                                          
%   * FillJf() => Jf_NZones_xxx() in MatLab                                                             
%     Diagonals:     dw/dPn                                                 
%     Off-Diagonals: dw/dPm                                                 
%                                                                           
%---------------------------------------------------------------------------
xeps = max(sim_params.xtol, 5*eps("double"));
feps = max(sim_params.ftol, 5*eps("double"));

if bitand(sim_params.verbose, logger.LOG_SIM_AF_DETAILS)
    nout = length(x0);
    fprintf('===== newtonSys.m =====\n');
    fprintf('xtol= %6.3e ftol= %6.3e relax= %4.2f maxit= %3d verbose= %1d\n\n', ...
        sim_params.xtol, sim_params.ftol, sim_params.relax, sim_params.maxit, sim_params.verbose);
    fprintf('eps = %6.3e xeps = %6.3e feps = %6.3e\n', eps("double"), xeps, feps);
end

%  Initial guess and current number of iterations
x = x0';
k = 0;        

Paths = varargin{1};
Zones = varargin{2};

% Calculate stack pressures for all airflow paths based on zone densities.
Paths = calc_SP( Paths, Zones, sim_params.verbose );

while k <= sim_params.maxit
    k = k + 1;

    % Fill the Jacobian matrix.
    %   Jfun, e.g., Jf_NZones_PowerLaw() returns the Jacobian matrix [J]
    %   and the residual vector {f}.
    [J,F,Paths,Zones] = feval( Jfun, sim_params, x, Paths, Zones);
    F = F';
    dx = J\F;
    x = x - sim_params.SYSTEM_SIGN * sim_params.relax * dx;

    if bitand(sim_params.verbose,logger.LOG_SIM_AF_DETAILS)
        fprintf('\n----- newtonSys -----\nNewton iterations\n  k     norm(f)      norm(dx) \n');
        fprintf('%3d %12.3e %12.3e \n', k, norm(F), norm(dx));
        for i=1:nout
            fprintf('    f(%1d)     ',i);
        end
        fprintf('\n');
        fprintf('%12.5e ', F);
        fprintf('\n');
        for i=1:nout
            fprintf('    x(%1d)     ',i);
        end
        fprintf('\n');
        fprintf('%12.5e ', x);
        fprintf('\n');
        for i=1:nout
            fprintf('     dx(%1d)   ',i);
        end
        fprintf('\n');
        fprintf('%12.5e ', dx);
        fprintf('\n');
    end
  
    varargout{1} = x;
    varargout{2} = Paths;
    varargout{3} = Zones;
    
    nF = norm(F);
    ndx = norm(dx);
  if ( nF < feps ) || ( ndx < xeps )
      if bitand(sim_params.verbose,logger.LOG_RESULTS)
          fprintf('\n*** Airflows Converged after %d iterations on ', k);
          if( nF < feps )
              fprintf( 'norm(f)= %.2e < feps= %.3e\n', nF, feps );
          else
              fprintf( 'norm(dx)= %.2e < xeps= %.3e\n', ndx, feps );
          end
      end
%      J
%      A = ones(size(J))
%      posdef = J'*A*J
      return;
  end
end
warning('Solution not found within tolerance after %d iterations\n',k);
end % end newtonSys()