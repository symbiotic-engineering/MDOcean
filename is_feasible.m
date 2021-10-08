function [feasible, failed] = is_feasible(power, B, FOS, p)

B_ok = B > p.B_min;
FOS_ok = FOS > p.FOS_min;
power_ok = power > p.P_min;

feasible = B_ok & FOS_ok & power_ok;

if nargout > 1
%     constraint_names = {'Buoyancy','FOS','Power'};
%     idx_failed = [];
%     if ~B_ok 
%         idx_failed = [idx_failed 1];
%     end
%     if ~FOS_ok
%         idx_failed = [idx_failed 2];
%     end
%     if ~power_ok
%         idx_failed = [idx_failed 3];
%     end
% 
%     which_constraints_failed = constraint_names(idx_failed);

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