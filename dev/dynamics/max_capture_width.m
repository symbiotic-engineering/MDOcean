syms T g
w = 2*pi/T; % definition of angular frequency
k = w^2/g; % dispersion relation for deep water
lambda = 2*pi/k; % definition of wavenumber
CW_max = lambda / (2*pi) % eq 3 https://www.sciencedirect.com/science/article/pii/S0960148115001652#bib50
