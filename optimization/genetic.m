clc, clear, close all
p = parameters();
b = var_bounds(p);

FitFcn=@(X)simulation(X,p);
nvars=7;
lb=[b.D_sft_min, b.D_i_ratio_min, b.D_or_ratio_min, b.M_min, b.N_WEC_min, b.D_int_min, b.w_n_min];
ub=[b.D_sft_max, b.D_i_ratio_max, b.D_or_ratio_max, b.M_max, b.N_WEC_max, b.D_int_max, b.w_n_max];
options = optimoptions('ga', 'PlotFcn', {@gaplotbestf}, 'Display', 'iter');%,'PopulationSize',100);
intcon = 4;

[X,fval,~,~,finalpop]=ga(FitFcn,nvars,[],[],[],[],lb,ub,@mycon,intcon,options);


%Constraints
function [c,ceq] = mycon(X)
p=parameters();
[~, ~, B, FOS1Y, FOS2Y, ...
            FOS3Y, FOS_buckling,GM, P_elec, ~] = simulation (X,p);

c(1)=-FOS1Y(1)+p.FOS_min;
c(2)=-FOS1Y(2)+p.FOS_min;%<=0;
c(3)=-FOS2Y(1)+p.FOS_min; %<=0;
c(4)=-FOS2Y(2)+p.FOS_min; %<=0;
c(5)=-FOS3Y(1)+p.FOS_min; %<=0;
c(6)=-FOS3Y(2)+p.FOS_min; %<=0;
c(7)=-FOS_buckling(1)+p.FOS_min; %<=0;
c(8)=-FOS_buckling(2)+p.FOS_min; %<=0;
c(9)=-B+p.B_min; %<=0;
c(10)= -GM;
c(11) = -P_elec;
ceq=[];
end