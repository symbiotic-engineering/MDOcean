function pts = conic_circle_intersect(a,b,c,d,e,f,h,k,R)
% INTERSECT A GENERAL CONIC WITH A CIRCLE
%
% Conic: a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0
% Circle: (x - h)^2 + (y - k)^2 = R^2
%
% Returns Nx2 array of intersection points [x, y]

%% Step 1: Rotate to remove xy term
theta = 0.5*atan2(b, a - c);
ct = cos(theta);
st = sin(theta);

% Rotation transform:
% x = ct*X - st*Y
% y = st*X + ct*Y

% Rotate conic coefficients
A_rot = a*ct^2 + b*ct*st + c*st^2;
C_rot = a*st^2 - b*ct*st + c*ct^2; % coefficient of Y^2
B_rot = 0; % xy term eliminated
D_rot = d*ct + e*st;
E_rot = -d*st + e*ct;
F_rot = f;

%% Step 2: Rotate circle center
H = h*ct + k*st;
K = -h*st + k*ct;

%% Step 3: Build quartic in X using Y from circle
% Circle: (X-H)^2 + (Y-K)^2 = R^2 -> Y = K +/- sqrt(R^2 - (X-H)^2)
%
% Conic: A_rot*X^2 + C_rot*Y^2 + D_rot*X + E_rot*Y + F_rot = 0
% Substitute Y

alpha2 = A_rot - C_rot;       % coefficient of X^2 inside square
alpha1 = D_rot + 2*C_rot*H;   % coefficient of X
alpha0 = C_rot*(K^2 + R^2 - H^2) + E_rot*K + F_rot; % constant term
beta0  = 2*C_rot*K + E_rot;   % multiplier of sqrt term

% Quartic: (alpha2*X^2 + alpha1*X + alpha0)^2 - beta0^2*(R^2 - (X-H)^2) = 0
p4 = alpha2^2;
p3 = 2*alpha2*alpha1;
p2 = alpha1^2 + 2*alpha2*alpha0 + beta0^2;
p1 = 2*alpha1*alpha0 - 2*beta0^2*H;
p0 = alpha0^2 - beta0^2*(R^2 - H^2);

quartic = [p4 p3 p2 p1 p0];

%% Step 4: Solve quartic
Xroots = roots(quartic);

% Keep real roots
Xroots = Xroots(abs(imag(Xroots)) < 1e-10);
Xroots = real(Xroots);

%% Step 5: Recover Y and rotate back
pts = [];
for i = 1:length(Xroots)
    X = Xroots(i);
    disc = R^2 - (X-H)^2;
    if disc < 0
        continue
    end
    Y_candidates = [K + sqrt(disc), K - sqrt(disc)];
    for Y = Y_candidates
        % Rotate back
        x =  ct*X - st*Y;
        y =  st*X + ct*Y;
        % Verify solution
        if abs(a*x^2 + b*x*y + c*y^2 + d*x + e*y + f) < 1e-8 && ...
           abs((x-h)^2 + (y-k)^2 - R^2) < 1e-8
            pts = [pts; x y];
        end
    end
end

% Remove duplicates
pts = unique(round(pts,12),'rows');

end
