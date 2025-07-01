function [m_dot, dFdP] = afe_plr(dP, path, verbose)
% Calculate mass flow etc for Powerlaw airflow elements
    
    % ContamX 3.4.x
    %pn->visc = (R4)(3.7143e-6 + 4.9286E-8 * pn->T);
    %pn->Dvisc = (R4)(pn->D / pn->visc);

    if(dP >= 0)
        T = path.Tn;
    else
        T = path.Tm;
    end
    T = T + 273.15;

    visc = 3.7143e-6 + 4.9286e-8 * T;   % kinematic [kg/m-s]
    Dvisc = path.D / visc;              % dynamic [s/m2]
    cdm = path.Clam * Dvisc;
    Fl = sign(dP) * cdm * abs(dP);
    Ft = sign(dP) * path.Cturb * sqrt(path.D) * realpow(abs(dP), path.expt);

    % Calculate derivative
    if( abs(Fl) <= abs(Ft) )
        m_dot = Fl;
        dFdP = cdm;
        Regime = 'lam';
    else
        m_dot = Ft;
        dFdP = path.expt * m_dot / dP;
        Regime = 'turb';
    end
    m_dot = path.mult * m_dot;
    dFdP = path.mult * dFdP;

    if bitand(verbose, logger.LOG_AFE)
        fprintf('  Path(%1d,%1d) viscK= %8.5e viscD= %8.5e kg/m-s dP= %+8.5e Pa   Flam= %+8.5e Fturb= %+8.5e mdot=%+8.5e kg/s regime=%s\n', ...
            path.zn, path.zm, visc, Dvisc, dP, Fl, Ft, m_dot, Regime );
    end

end