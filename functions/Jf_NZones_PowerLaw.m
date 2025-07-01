function [J, f, Paths, Zones] = Jf_NZones_PowerLaw(sim_params, x, Paths, Zones)
% Jf_NZones_PowerLaw()  Jacobian and f vector for generic n-zones system
%
% Synopsis:  [J, f, Paths, Zones] = Jf_NZones_PowerLaw(sim_params, x, Paths, Zones)
%
% Input: sim_params = simulation parameters 
%        x          = current guess at solution vector
%        Paths      = vector of Paths structures
%        Zones      = vector of Zones structures
%
% Output: J     = Jacobian matrix for the non-linear system
%         f     = right hand side vector for Newton iterations
%         Paths = vector of airflow path structures
%         Zones = vector of zone structures, density will change based on pressures

DEBUG = 2;

nz = length(x);         % Number of pressure nodes (or zones not including ambt) 
np = length(Paths);     % Number of airflow paths

% Extract problem variables, i.e., zone pressures, from x,
%   and set ambient pressure to zero.
p  = [0;x];             
% Reset residuals and Jacobian to zero.
f  = zeros(1,nz);       % f = residual mass flow equations
J  = zeros(nz,nz);      % J = Jacobian matrix

% Calculate zone densities.
for iz=1:nz+1
    Zones(iz).D = Density(101325 + p(iz) * sim_params.USE_ZONE_PRESSURE, Zones(iz).T);
    if(DEBUG > 3)
        fprintf("Zone(%d): P = %.5e RHO = %7.5f\n", iz, p(iz), Zones(iz).D)
    end
end
% Calculate stack pressures.
Paths = calc_SP( Paths, Zones, sim_params.verbose );

if bitand(sim_params.verbose, logger.LOG_SIM_AF)
    fprintf('path  zm  zn     pm         pn           D         dP            F          dFdP\n');
end
for i=1:np
    % Determine zones m and n for path i
    m = Paths(i).zm;
    n = Paths(i).zn;

    % Calculate dP = pm - pn + Pwind + Pstack
    pm = p(m);
    pn = p(n);
    dP = pn - pm + Paths(i).Pw + Paths(i).Ps;    
    
    % Set path density to that of incoming air OR average.
    if( strcmp(sim_params.path_density,'AVERAGE') )
        rho = 0.5*(Zones(m).D + Zones(n).D);
    else
        if( dP > 0.0 )
            rho = Zones(n).D;
        else
            rho = Zones(m).D;
        end
    end
    Paths(i).D = rho;
    
    % Calculate mass flow rates and derivatives for all flow paths.
    % Currently, only powerlaw elements are implemented, and all are using
    % the same element type.
    [F, dFdP] = afe_plr(dP, Paths(i), sim_params.verbose);
    
    % Fill Jacobian
    % Diagonal
    %   skip paths connected to constant pressure zones, e.g., 1
    if( n ~= 1 )
        J(n-1,n-1) = J(n-1,n-1) + dFdP;
        f(Paths(i).zn-1) = f(Paths(i).zn-1) + sim_params.SYSTEM_SIGN * F;
    end
    if( m ~= 1 )
        J(m-1,m-1) = J(m-1,m-1) + dFdP;
        f(Paths(i).zm-1) = f(Paths(i).zm-1) - sim_params.SYSTEM_SIGN * F;
    end

    % Off-diagonals
    %   skip paths connected to constant pressure zones, e.g., zone 1
    if( m ~=1 && n ~= 1 )
        J(m-1,n-1) = J(n-1,m-1) - dFdP;
        J(n-1,m-1) = J(n-1,m-1) - dFdP;
    end

    if bitand(sim_params.verbose, logger.LOG_SIM_AF)
        % Log iterate.
        fprintf('%4d %3d %3d %+8.3e %+8.3e  %8.4f  %+8.5e  %+8.5e %+8.5e\n', ...
            i, m, n, pm, pn, rho, dP, F, dFdP);
    end
    
end

if bitand(sim_params.verbose, logger.LOG_SIM_AF_DETAILS)
    J
    F
    Zones
end
