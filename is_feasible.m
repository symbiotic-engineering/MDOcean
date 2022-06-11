function [feasible, failed] = is_feasible(B, FOS, GM, P_elec, D_d, h_s_extra, LCOE_const, F_max_const, p)

B_ok = all(B < p.B_min) & all(B > 0);
FOS_ok = FOS > p.FOS_min;
GM_ok = GM > 0;
P_ok = P_elec > 0;
D_d_ok = D_d >= p.D_d_min;
h_s_ex_ok = h_s_extra >= 0;
L_ok = LCOE_const >= 0;
F_ok = F_max_const >= 0;

consts = [B_ok FOS_ok GM_ok P_ok D_d_ok h_s_ex_ok L_ok F_ok];
feasible = all(consts);
const_names = {'Buoyancy ','FOS ','Metacentric Height ','Power ','Damping ','Spar Height ','LCOE max ','F max '};


if nargout > 1
    failed = ' ';
    for i=1:length(consts)
        if ~consts(i)
            failed = [failed const_names{i}];
        end
    end
end

end