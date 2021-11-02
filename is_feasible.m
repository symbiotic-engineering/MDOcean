function [feasible, failed] = is_feasible(B, FOS, GM, p)

B_ok = B > p.B_min;
FOS_ok = FOS > p.FOS_min;
GM_ok = GM > 0;

feasible = B_ok & FOS_ok & GM_ok;

if nargout > 1
    failed = ' ';
    if ~B_ok 
        failed = [failed 'Buoyancy '];
    end
    if ~FOS_ok
        failed = [failed 'FOS '];
    end
    if ~GM_ok
        failed= [failed 'Metacentric Height'];
    end
end

end