function [feasible, failed] = is_feasible(power, B, FOS, p)

B_ok = B > p.B_min;
FOS_ok = FOS > p.FOS_min;
power_ok = power > p.P_min;

feasible = B_ok & FOS_ok & power_ok;

if nargout > 1
    failed = ' ';
    if ~B_ok 
        failed = [failed 'Buoyancy '];
    end
    if ~FOS_ok
        failed = [failed 'FOS '];
    end
    if ~power_ok
        failed = [failed 'Power'];
    end  
end

end