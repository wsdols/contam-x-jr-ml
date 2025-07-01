function [Cturb, Clam] = setFlowCoef(Path, element_type)
%---------------------------------------------------%
% Calculate Turbulent and Laminar flow coefficients %
%---------------------------------------------------%
% Re = rho * V * D / mu;
% F = rho * V * A

MUAIR  = 1.81625e-05;       % Dynamic viscosity
DPTMIN = 1.0e-10;
RHOAIR = 1.20410;
SQRTRHO = sqrt(RHOAIR);
SQRT2 = sqrt(2.0);
ReTrans = 30.0;             % For PL_LEAK elements

if(strcmp(element_type,'PL_ORFC'))
    % F = Cd * A sqrt(2*rho*dP)
    ReTrans = Path.ReTr;
    Cturb = Path.Cd * Path.area * SQRT2;
    Ftrans = MUAIR * ReTrans * Path.area / sqrt(4 * Path.area / 3.14159);
elseif(strcmp(element_type,'PL_LEAK3'))
    dia = sqrt(Path.area);
    Cturb = SQRT2 * Path.Cd * Path.area * power(Path.Pref, 0.5 - Path.expt);
    Ftrans = MUAIR * ReTrans * dia;
else
    fprintf('ERROR - %s is not recognized\n', element_type);
end
dPtrans = Ftrans / (Cturb*SQRTRHO);
dPtrans = power(dPtrans,1.0/Path.expt);

if( dPtrans < DPTMIN ) 
  dPtrans = DPTMIN; 
end
Clam = (MUAIR * Ftrans) / (RHOAIR * dPtrans);

% F = Cl * (rho / mu) * dP

return