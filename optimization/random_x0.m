function [x0,x0_struct] = random_x0(b)

mins = b.X_mins;
maxs = b.X_maxs;

r = rand(1,7);

x0 = (maxs-mins).*r + mins;
x0(8) = randi(b.M_max);

x0_struct = struct('D_f',x0(1),'D_s_ratio',x0(2),'h_f_ratio',x0(3),...
    'T_s_ratio',x0(4),'F_max',x0(5),'D_int',x0(6),'w_n',x0(7));

end

