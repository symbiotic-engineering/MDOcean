function trimmed = trim_jpd(jpd)

jpd_entries = jpd(2:end,2:end);
rows_all_zeros = all(jpd_entries==0,2);
cols_all_zeros = all(jpd_entries==0,1);
rows_keep = [true; ~rows_all_zeros];
cols_keep = [true ~cols_all_zeros];
trimmed = jpd(rows_keep, cols_keep);

min_T = trimmed(1,2);
if min_T < 3.5
    error(['The series expansion used to compute the froude krylov force' ...
        'coefficient in dynamics.m requires 2 / (T * R^(5/8)) << 1.' ...
        'If pushing up against this limit, it is recommended to plot ' ...
        'gamma vs omega to confirm smooth behavior at high frequencies. '])
end

end