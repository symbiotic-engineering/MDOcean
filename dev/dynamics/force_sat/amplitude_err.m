% notebook p113 2/7/25
% this is the amplitude ratio assuming zeta and wn vs freq are const
% I would need to combine this with the wave spectrum to get % energy accounted for
% and technically also combine with zeta and wn vs freq and excitation force vs freq

syms a positive;
syms b real;
syms n integer positive;
trm = 1 / sqrt(a * (1 + 2 * n)^4 + b * (1 + 2 * n)^2 + 1);
s = symsum(trm, n, 1, Inf);
simplify(s, 'Steps', 100);
