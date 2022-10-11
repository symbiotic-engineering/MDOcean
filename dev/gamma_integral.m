syms w t k rho R a T g real positive
syms x y real
phi_I = a*w/k * exp(-k*T) * 1i * exp(1i*(w*t - k*x));
integrand_pressure = rho * diff(phi_I,t);
integrand_gamma = simplify(integrand_pressure / (a*exp(1i*w*t)));

x_min = -sqrt(R^2-y^2);
x_max = sqrt(R^2-y^2);
inner_int = int(integrand_gamma,x,x_min,x_max);

coeffs = rho*w^2*exp(-T*k)/k^2;
inner_int_without_coeffs = inner_int / coeffs;

y_min = -R;
y_max = R;
gamma_without_coeffs = int(inner_int_without_coeffs,y,y_min,y_max);
gamma_without_coeffs = simplify(rewrite(gamma_without_coeffs,'sin'));
coeffs = subs(coeffs,k^2,k*w^2/g);
gamma = abs(gamma_without_coeffs * coeffs)
taylor(gamma_without_coeffs/k,k,'Order',11)

derivative = diff(gamma_without_coeffs,'k')

%% try hilbert transform
clear all
syms rho g k T R w real positive
bessel_approx = taylor(besselj(1,k*R),k,'Order',7);
gamma_bessel = 2*rho*g/k * exp(-k*T) * R * pi * bessel_approx;
b_over_rho_w = k/2 * (gamma_bessel/(rho*g))^2;
b = b_over_rho_w * rho * w;

b = subs(b,k,w^2/g);
pretty(b)


a = 1/w^2 * htrans(-w*b,w,w);
pretty(a)

fudge = 5;
b_over_rho_w_eval = vpa(subs(b_over_rho_w,{R,T,k},{10,2*fudge,w^2/9.8}));
a_over_rho_eval = vpa(subs(a/rho,{g,T,R,rho},{9.8,2*fudge,10,1000}));
pretty(a_over_rho_eval)

figure
fplot(a_over_rho_eval,[.1 2])
title('a/rho')
figure
fplot(b_over_rho_w_eval,[.1 1.5])
title('b/rho w')

%% approximation error
i = 8:2:14; % last term included
exc = i+2;  % first term excluded
w = 1.4;
k = w^2/9.8;
R_max = 20;
order_err = (k*R_max).^exc ./ factorial(exc);
figure
plot(i,order_err*100)
xlabel('exponent of k on last term included')
ylabel('Approx percent error, 0-100% scale')

%%


% syms r theta
% integrand_gamma_polar = r * subs(integrand_gamma,x,r*cos(theta));
% inner_int = int(integrand_gamma_polar,r,0,R);
% gamma = int(inner_int, theta, 0, 2*pi)

% subscript _num is for numerical, so I don't overwrite symbolic vars
R_num = [3 6 10 12 15 18];
T_wave = 5:2:13;
w_num = 2*pi./T_wave;
[R_mesh,w_mesh] = meshgrid(R_num,w_num);
k_mesh = w_mesh.^2/9.8;
approx = 2.88*R_mesh.^2;
actual = 1./k_mesh .* vpa(subs(gamma_without_coeffs,{R,k},{R_mesh,k_mesh}))

%%
w_mesh = sqrt(k_mesh*9.8);
actual=abs(actual);
error = actual./approx;

figure
subplot 121
contourf(R_mesh,k_mesh,approx)
colorbar
title('approx')
hold on
plot([10 10],[min(k_num) max(k_num)],'w--')

subplot 122
contourf(R_mesh,k_mesh,actual)
colorbar
title('actual')
hold on
plot([10 10],[min(k_num) max(k_num)],'w--')

figure
subplot 121
contourf(R_mesh,k_mesh,error)
colorbar
caxis([0 1])
title('error ratio')
hold on
plot([10 10],[min(k_num) max(k_num)],'w--')

subplot 122
trial = 1./(k_mesh .* R_mesh); %exp(-w_mesh.^2*0.5
contourf(R_mesh,k_mesh,trial/2)
colorbar
caxis([0 1])
%%
% line for each R
figure
subplot 121
for i=1:length(R_num)
plot(w_num,error(:,i),'DisplayName',num2str(R_num(i)))
hold on
end
title('Error for various radii')
xlabel('k')
legend

% line for each k
subplot 122
for i=1:length(k_num)
plot(R_num,error(i,:),'DisplayName',num2str(k_num(i)))
hold on
end
title('Error for various k')
xlabel('R')
legend

%%
% compare numerical integral to actual wamit
w_wamit = hydro.w;
gamma = hydro.ex_ma(3,1,:);
gamma = gamma(:)';
draft = 2;
wamit_compare = gamma ./ exp(-w_wamit.^2*draft/9.8)
mdocean = 314 * exp(w_wamit.^2*draft/9.8);

figure
plot(w_num,actual(:,R_num==10)-100)
hold on
plot(w_num,approx(:,R_num==10))
plot(w_wamit,mdocean)
plot(w_wamit,wamit_compare)
legend('Numerical integral','Approximate Integral (Desmos)','MDOcean original','WAMIT')
xlim([min(w_num) max(w_num)])
xlabel('w')
ylabel('$\frac{\gamma}{\rho g e^{-kT}}$','Interpreter','latex','FontSize',20)

%%
figure
plot(w_wamit,gamma)
hold on
fudge_factor = 1;
plot(w_num,actual(:,R_num==10) .* exp(-k_num'*fudge_factor*draft))
legend('WAMIT','Numerical Integration with 5x fudge')
xlim([min(w_num) max(w_num)])
xlabel('w')
ylabel('\gamma/\rho g')
