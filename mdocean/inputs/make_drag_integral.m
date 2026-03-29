drag_fcn = make_drag_integral_LUT(.01, .5, 50.01, pi/24, .1, 8);
save('./mdocean/inputs/drag_integral.mat','drag_fcn');
