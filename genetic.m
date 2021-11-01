clc, clear, close all
p = parameters();
p.M = 1; % hardcode steel for now - run a for-loop

FitFcn=@(X)simulation([X(1) X(2) X(3) 0 X(4) X(5)],zeros(20,21),p);
nvars=5;
lb=[0,0,0,1,1e6];
ub=[50,1,50,100,1e8];
options = optimoptions('ga', 'PlotFcn', {@gaplotbestf}, 'Display', 'iter');
[X,fval]=ga(FitFcn,nvars,[],[],[],[],lb,ub,NONLCON,options);

%Constraints
FOS1Y(1) 


