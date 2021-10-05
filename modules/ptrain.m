function [F_ptrain,P_elec] = ptrain(s, u, i_PT)

GR = 10; % gear ratio
r = 0.1;

T_max = 100; % max motor torque (in reality this is a function of i_PT)
F_max = T_max / r * GR;
F_ptrain = min(u,F_max);

T = F_ptrain * r / GR;      % motor torque Nm
w = s(2,:)/r * GR;          % motor speed rad/s
w_rpm = w*60/(2*pi);        % motor speed rpm

a = 2;
b = .00003;
c = 0;
d = 300;

P_mech = T.*w;
Loss = 0;%a*T.^2 + b*w_rpm.^3 + c*w_rpm + d;
P_elec = P_mech - Loss;

end

