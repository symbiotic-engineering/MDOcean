function [feasible, failed] = is_feasible(B, FOS, GM, P_elec, p)

B_ok = all(B < p.B_min) & all(B > 0);
FOS_ok = FOS > p.FOS_min;
GM_ok = GM > 0;
P_ok = P_elec > 0;

feasible = B_ok & FOS_ok & GM_ok & P_ok;

if nargout > 1
    failed = ' ';
    if ~B_ok 
        failed = [failed 'Buoyancy '];
    end
    if ~FOS_ok
        failed = [failed 'FOS '];
    end
    if ~GM_ok
        failed= [failed 'Metacentric Height '];
    end
    if ~P_ok
        failed = [failed 'Power'];
    end
end

end