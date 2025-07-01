function [Pz, mdot, dP] = airnet_pl2(log_level, Mult)
arguments
    % Set log_level to 1 if running directly.
    log_level(1,1) double {mustBeInteger} = 1;
    % Set Mult to 1:1,000,000:10 for 7 cases.     
    Mult(1,1) double {mustBeFloat, mustBePositive} = 1;
end
% airnet_pl1
% Version: 6/12/2025
%
% Synopsis:
%   AIRNET Powerlaw Element Test # 1
%
% Output:  Solution is printed to command window

verbose = log_level;

% ---------- Simulation parameters ----------
sim_params = struct( ...
    'verbose', verbose,...  % See logger.m for BIT MASKS
    'xtol', 1e-9,...        % Simulation tolerance on pressure vector {x}
    'ftol', 1e-9,...        % Simulation tolerance on residual vector {F} 
    'relax', 0.75,...       % Relaxation factor (GNW => 0.75 works well)
    'maxit', 30,...         % Maximum number of non-linear iterations
    ... SYSTEM_SIGN - See newtonSys() comments:
    ...   +1 => As done in ContamX
    ...   -1 => As formulated in Dennis and Schnabel
    'SYSTEM_SIGN', -1,...
    ... Path density to use: 'DIRECTIONAL' or 'AVERAGE'
    'path_density','DIRECTIONAL',...
    ... Zone density calculation method:
    ...   0 => use barametric pressure, 
    ...   1 => use zone pressure
    'USE_ZONE_PRESSURE',0,...
    ... Airflow element type to use for airflow paths (all the same)
    ...    'PL_ORFC', 'PL_LEAK3'
    'element_type', 'PL_ORFC'...
    );

% ---------- Weather, i.e., ambient parameters ----------
wth = struct( ...
    'Pambt',101325.0,...	% Pa
    'Tambt',20.0,...        % C
    'WndSpd',0.0,...        % m/s
    'Rair',287.055...       % J/kg-K
    );
Dambt = wth.Pambt/(wth.Rair * (wth.Tambt+273.15));  % Ambt density

%==========================================================================
% Project Configuration
%==========================================================================
% This model represents test case found in the AIRNET NISTIR 89-4072:
%    Powerlaw Element Test #2.
%    Results are slightly different due to changes in CONTAM.
%    Eight separate zones (A-F) with different wind pressure boundaries.
%    "Very" low pressures result in laminar flows for cases F-H. 
%
% CONTAM z#-name
%      ┌──────────┬──────────┐ 
%  1   ◊          │          │ 
%      │          ◊          │ 
% ambt │ 2,A      │ 3,B      ◊ 
%      └──────────┴──────────┘ 
%==========================================================================
% Create array of Zone structures
% NUMBER (CONTAM #)       1(0)     2(1) 3(2) 
Zones = struct('name', { 'ambt'    'A'  'B'  },...
               'T',    { wth.Tambt 20.0 20.0 },...
               'Zref', {   0.0      0.0  0.0 },...
               'P0',   {   0.0      0.0  0.0 },...
               'D',    {   0.0      0.0  0.0 }...
              );
Aorf1 = 1e-4;   % Orifice area [m2]
Pw1= 100.;      % Wind pressure [Pa]

% Create array of Path structures
% CONTAM #    1     2     3    
Paths = struct( ... 
    'zn',   { 1     2     1    },...      % Positive flow from n to m
    'zm',   { 2     3     3    },...
    'area', { Aorf1 Aorf1 Aorf1},...      % Opening area m2/m2
    'Cd',   { 0.6   0.6   0.6  },...      % Discharge coefficients
    'expt' ,{ 0.5   0.5   0.5  },...      % Use with PL_LEAK
    'Pref' ,{ 4.0   4.0   4.0  },...      % Use with PL_LEAK
    'ReTr' ,{100.  100.  100.  },...      % Transition Reynolds number
    'mult', { 1.0  Mult   1.0  },...
    'Clam', { 0.0   0.0   0.0  },...      % Laminar flow coefficient
    'Cturb',{ 0.0   0.0   0.0  },...      % Turbulent flow coefficient
    'Z',    { 0.0   0.0   0.0  },...      % Absolute height of path (m)
    'Ht_m', { 0.0   0.0   0.0  },...
    'Ht_n', { 0.0   0.0   0.0  },...
    'Tm',   { 0.0   0.0   0.0  },...
    'Tn',   { 0.0   0.0   0.0  },...
    'Cp',   { 1.0   1.0   1.0  },...
    'Ps',   { 0.0   0.0   0.0  },...      % Initialize stack pressures (Pa)
    'Pw',   { Pw1   0.0   0.0  },...
    'D',    { 0.0   0.0   0.0  }...
    );
%==========================================================================
              
% finish setting Zone and Path data
nz = length(Zones);
for i=1:nz
    Zones(i).D = Density( wth.Pambt, Zones(i).T );
end
if bitand(sim_params.verbose, logger.LOG_PROJECT)
fprintf('\n=================================================== INPUTS =====\n\n');
    fprintf('WEATHER DATA\n')
    fprintf('\tPambt:\t%.0f Pa\n\tTambt:\t%6.2f C\n\tWind Speed:\t%4.2f m/s\n', ...
        wth.Pambt, wth.Tambt, wth.WndSpd)
 
    fprintf('ZONE DATA\n')
    fprintf('  #   Name     T    Rho(initial)\n')
    for i=1:nz
        fprintf('%3d %6s %8.3f %7.5f\n', i, Zones(i).name, Zones(i).T, Zones(i).D)
    end
end

if bitand(sim_params.verbose, logger.LOG_PROJECT)
    fprintf('PATH DATA\n');
    fprintf(' #   zm  zn   Tm     Tn   area[m2]   mult        Cp      Pw          Cturb        Clam        Ht_m   Ht_n\n');
end
np = length(Paths);
for i=1:np
    % Set Path inlet and outlet temperatures
    Paths(i).Tn = Zones(Paths(i).zn).T;
    Paths(i).Tm = Zones(Paths(i).zm).T;
    % Calculate Path wind pressures (Pa)
    if(abs(Paths(i).Pw) == 0.0)
        Pw = 0.5 * Paths(i).Cp * Dambt * wth.WndSpd * wth.WndSpd;
        Paths(i).Pw = Pw;
    end
    
    % Calculate turbulent and laminar flow coefficients of flow path
    [Paths(i).Cturb, Paths(i).Clam] = setFlowCoef(Paths(i), sim_params.element_type);

    % Calculate absolute height and Initial density for each path
    Paths(i).Ht_m = Paths(i).Z - Zones(Paths(i).zm).Zref;
    Paths(i).Ht_n = Paths(i).Z - Zones(Paths(i).zn).Zref;
    Paths(i).D = 0.5 * ( Density(wth.Pambt,Paths(i).Tm) + ...
                         Density(wth.Pambt,Paths(i).Tn) );

    if bitand(sim_params.verbose, logger.LOG_PROJECT)
        fprintf(['%3d %3d %3d %6.2f %6.2f ' ...
            '%8g %7.3e %8.3f %11.4e  ' ...
            '%11.4e  %11.4e %6.2f %6.2f\n'], ...
            i, Paths(i).zm, Paths(i).zn, Paths(i).Tm, Paths(i).Tn, ...
            Paths(i).area, Paths(i).mult, Paths(i).Cp, Paths(i).Pw, ...
            Paths(i).Cturb, Paths(i).Clam, Paths(i).Ht_m, Paths(i).Ht_n);
    end
end

if bitand(sim_params.verbose, logger.LOG_RESULTS)
    fprintf('\n================================================== RESULTS =====\n');
end

% --- Solve the system and print results
P0   = zeros(1,nz-1);   % n-1 vector of unknowns
x0 = P0;                %  initial guess at solution vector
[x, Paths, Zones] = newtonSys('Jf_NZones_PowerLaw', sim_params, x0, Paths, Zones);
Pz = [0;x];
[mdot, dP] = CalculateAllFlows(Pz, Paths);

if bitand(sim_params.verbose, logger.LOG_RESULTS)
    fprintf('\nNOTE: Diagram appears best with Courier New font.\n\n')
    fprintf('      ┌──────────┬──────────┐  \n');
    fprintf('  1   ◊          │          │  \n');
    fprintf('      │          ◊          │  \n');
    fprintf(' ambt │ 2,A      │ 3,B      ◊  \n');
    fprintf('      └──────────┴──────────┘  \n');
    fprintf('\n')
    
    fprintf('Tambt %5.2f  WindSpeed %5.2f\n', wth.Tambt, wth.WndSpd)
    fprintf('"Zone Results"\n');
    fprintf('  #   Name   Pressure[Pa]  Density[kg/m3]\n')
    for i=1:nz
        fprintf('%3d %6s   %+9.4e   %.5f\n', i, Zones(i).name, Pz(i), Zones(i).D );
    end
    fprintf('\n"Airflow Path Results"\n');
    fprintf('  #    mdot[kg/s]   dP[Pa]      D[kg/m3]   Pw[Pa]      Ps[Pa]     FlowDir\n')
    
    for i=1:np
        if(dP(i) > 0)
            % From Zone n To Zone m
            flow_dir_str = num2str(Paths(i).zm) + "<--" + num2str(Paths(i).zn);
        else
            % From Zone m To Zone n
            flow_dir_str = num2str(Paths(i).zm) + "-->" + num2str(Paths(i).zn);
        end
        fprintf('%3d  %+8.5e %+8.5e %8.5f  %+9.3e  %+9.3e    %s\n', ...
            i, mdot(i), dP(i), Paths(i).D, Paths(i).Pw, Paths(i).Ps, flow_dir_str);
    end
end

end % end airnet_pl2()