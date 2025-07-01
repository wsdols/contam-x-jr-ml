function rho = Density(P,T)

Rair  = 287.055;                    % Universal gas constant
rho   = P/(Rair * (T+273.15));      % Density
