function [x0,x0_struct] = random_x0(b)

mins = b.X_mins;
maxs = b.X_maxs;

r = rand(length(mins),1);

x0 = (maxs-mins).*r + mins;
x0(length(b.var_names)) = randi(b.M_max);

x0_struct = cell2struct(num2cell(x0(1:end-1)),b.var_names(1:end-1));

end
