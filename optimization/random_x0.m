function [x0,x0_struct] = random_x0(b)

mins = [b.D_sft_min, b.D_i_ratio_min, b.D_or_ratio_min, b.M_min, b.N_WEC_min, b.D_int_min, b.w_n_min];
maxs = [b.D_sft_max, b.D_i_ratio_max, b.D_or_ratio_max, b.M_max, b.N_WEC_max, b.D_int_max, b.w_n_max];

r = rand(1,7);

x0 = (maxs-mins).*r + mins;
x0(4) = randi(b.M_max);

x0_struct = struct('D_sft',x0(1),'D_i_ratio',x0(2),'D_or_ratio',...
    x0(3),'N_WEC',x0(5),'D_int',x0(6),'w_n',x0(7));

end

