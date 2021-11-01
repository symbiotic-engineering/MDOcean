clc, clear, close all
p = parameters();
p.M = 1; % hardcode steel for now - run a for-loop

FitFcn=@(X)simulation([X(1) X(2) X(3) 0 X(4) X(5)],zeros(20,21),p);
nvars=5;
lb=[0,0,0,1,1e6];
ub=[50,1,50,100,1e8];
options = optimoptions('ga', 'PlotFcn', {@gaplotbestf}, 'Display', 'iter');
[X,fval]=ga(FitFcn,nvars,[],[],[],[],lb,ub,@mycon,options);


%Constraints
function [c,ceq] = mycon(X)
p=parameters();
[LCOE, D_env, B, FOS1Y, FOS2Y, ...
            FOS3Y, FOS_buckling, ...
            ~, ~, F_pt_unsat] = simulation ([X(1) X(2) X(3) 0 X(4) X(5)],zeros(20,21),p);

c(1)=FOS1Y(1)-p.FOS_min; %<=0;
c(2)=FOS1Y(2)-p.FOS_min;%<=0;
c(3)=FOS2Y(1)-p.FOS_min; %<=0;
c(4)=FOS2Y(2)-p.FOS_min; %<=0;
c(5)=FOS3Y(1)-p.FOS_min; %<=0;
c(6)=FOS3Y(2)-p.FOS_min; %<=0;
c(7)=FOS_buckling(1)-p.FOS_min; %<=0;
c(8)=FOS_buckling(2)-p.FOS_min; %<=0;
c(9)=B-p.B_min; %<=0; 
ceq=[];
end