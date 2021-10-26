function [feasible, failed] = is_feasible(power, B, FOS, GM, p)

B_ok = B > p.B_min;
FOS_ok = FOS > p.FOS_min;
power_ok = power > p.P_min;
GM_ok=GM > 0;

feasible = B_ok & FOS_ok & power_ok & GM_ok;

if nargout > 1
%     constraint_names = {'Buoyancy','FOS','Power','Metacentric Height'};
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
%     if ~GM_ok
%         idx_failed = [idx_failed 4};
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
    if ~GM_ok
        failed= [failed 'Metacentric Height'];
    end

end

end