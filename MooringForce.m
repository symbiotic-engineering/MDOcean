% Mooring Force Calculations
clc
syms T1 T2 F_moor theta2 theta1
% L_1 = 27.25; % length of the chain
% L_2 = 231.4; % length of nylon rope to sea buoy
theta1 = 20;
multiplier = 10;
B = multiplier*55000;
W = 98394.25;

sumFx_b = T2*cosd(theta2)==F_moor;
sumFy_b = T2*sind(theta2)==B;
sumFy_c = T1*sind(theta1)+W==T2*sind(theta2);
sumFx_c = T1*cosd(theta1)==T2*cosd(theta2);
% sumMz_a = W*(L_1*cosd(theta1))+F_moor*(L_1*sind(theta1)+L_2*sind(theta2))==B*(L_1*cosd(theta1)+L_2*cosd(theta2));

S=solve(sumFx_b,sumFy_b,sumFy_c,sumFx_c,T1,T2,F_moor,theta2);
F = eval(S.F_moor)
th2 = eval(S.theta2)

T1 = eval(S.T1)
T2 = eval(S.T2)