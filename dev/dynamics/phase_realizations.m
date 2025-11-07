% page 43 notebook WecSim results
% see slide 66 and 94 in summer 24 slideshow
powers = [13116 8658 13346 14054 13441 12216 13135 13839 15138 10817];

means = zeros(size(powers));
for i=1:length(powers)
    power_so_far = powers(1:i);
    means(i) = mean(power_so_far);
end

figure
plot(1:length(powers),means)
xlabel('Number of phase realizations')
ylabel('Mean power')
