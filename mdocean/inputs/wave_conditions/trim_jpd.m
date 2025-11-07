function trimmed = trim_jpd(jpd)

jpd_entries = jpd(2:end,2:end);
rows_all_zeros = all(jpd_entries==0,2);
cols_all_zeros = all(jpd_entries==0,1);
rows_keep = [true; ~rows_all_zeros];
cols_keep = [true ~cols_all_zeros];
trimmed = jpd(rows_keep, cols_keep);

end
