function c_v = nominal_c_v()
	filename = 'RM3-CBS.xlsx'; % spreadsheet containing RM3 "actual" power data
	sheet = 'Performance & Economics'; % name of relevant sheet
	
	P_matrix = readmatrix(filename,'Range','E97:S110','Sheet',sheet);
	JPD = readmatrix(filename,'Range','E24:S37','Sheet',sheet);
	
	P_weighted = P_matrix .* JPD / 100;
	P_elec = sum(P_weighted(:));

	% coefficient of variance (normalized standard deviation) of power
	P_var = std(P_matrix(:), in.JPD(:)) / P_elec;
	P_var = P_var * 100; % convert to percentage

end