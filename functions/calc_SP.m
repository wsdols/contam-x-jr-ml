function Paths = calc_SP( Paths, Zones, verbose )
% function npath = calc_SP()

if bitand(verbose, logger.LOG_SIM_AF)
    fprintf('calc_SP()\n')
    fprintf('path       Zn         Z1         h1         Dn         Zm         Z2         h2         Dm         dPs\n' );
end

g = 9.8055;
np = length(Paths);     % Number of airflow paths
for i=1:np
    % Determine zones m and n for path i
    m = Paths(i).zm;
    n = Paths(i).zn;
    Zm = Zones(m).Zref;
    Zn = Zones(n).Zref;
    h1 = Paths(i).Ht_m;
    h2 = Paths(i).Ht_n;
    Dm = Zones(m).D;
    Dn = Zones(n).D;

    % based on average density
    Ps = g * (Dm*h1 - Dn*h2);
    Paths(i).Ps = Ps;
    
    if bitand(verbose, logger.LOG_SIM_AF)
        fprintf( '%3d%c %10.4f %10.4f %10.4f %10.5f %10.4f %10.4f %10.4f %10.5f %10.6f\n',...
                 i, 'p', Zn, Paths(i).Z, h2, Dn, Zm, Paths(i).Z, h1, Dm, Paths(i).Ps);
    end
end
end % end calc_SP()
