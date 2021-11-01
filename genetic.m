clc, clear, close all
p = parameters();

w_max = 2*pi/p.T(find(any(p.JPD > 0),1,'first')); % max wave frequency that has any energy
w_min = 2*pi/p.T(find(any(p.JPD > 0),1,'last'));  % min wave frequency that has any energy

FitFcn=@(X)simulation(X,p);
nvars=7;
lb=[0,0,0,1,1,1,w_min];
ub=[50,1,50,3,100,100,w_max];
options = optimoptions('ga', 'PlotFcn', {@gaplotbestf}, 'Display', 'iter');
intcon = 4;
[X,fval]=ga(FitFcn,nvars,[],[],[],[],lb,ub,@mycon,intcon,options);


%Constraints
function [c,ceq] = mycon(X)
p=parameters();
[LCOE, D_env, B, FOS1Y, FOS2Y, ...
            FOS3Y, FOS_buckling,GM ...
            ~, ~, F_pt_unsat] = simulation (X,p);

c(1)=FOS1Y(1)-p.FOS_min; %<=0;
c(2)=FOS1Y(2)-p.FOS_min;%<=0;
c(3)=FOS2Y(1)-p.FOS_min; %<=0;
c(4)=FOS2Y(2)-p.FOS_min; %<=0;
c(5)=FOS3Y(1)-p.FOS_min; %<=0;
c(6)=FOS3Y(2)-p.FOS_min; %<=0;
c(7)=FOS_buckling(1)-p.FOS_min; %<=0;
c(8)=FOS_buckling(2)-p.FOS_min; %<=0;
c(9)=B-p.B_min; %<=0;
c(10)=GM;
ceq=[];
end