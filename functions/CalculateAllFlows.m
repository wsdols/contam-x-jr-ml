function [m_dot, dP] = CalculateAllFlows(p, paths)
% Calculate mass flow rate and pressure drop across flow path for
%   powerlaw airflow elements.
%
% Input:  
%   p     = current zone pressures
%   paths = array of Path structures
% Output: 
%   m_dot = mass flow rate through each Path
%   dP    = pressure difference accros each Path
%           +'ve defined as from n to m

N = length(paths);
m_dot = zeros(1,N);
dP    = zeros(1,N);
for i=1:N
    dP(i)    = p(paths(i).zn) - p(paths(i).zm) + paths(i).Pw + paths(i).Ps;
    [m_dot(i),~] = afe_plr(dP(i),paths(i),0);
%    m_dot(i) = sign(dP(i)) * paths(i).mult * paths(i).Cturb * sqrt(paths(i).D)...
%               * power(abs(dP(i)), paths(i).expt);
end