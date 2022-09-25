function c_v = nominal_c_v()
	filename = 'RM3-CBS.xlsx'; % spreadsheet containing RM3 "actual" power data
	sheet = 'Performance & Economics'; % name of relevant sheet
	
	P_matrix = readmatrix(filename,'Range','E97:S110','Sheet',sheet) * 1000;
	JPD = readmatrix(filename,'Range','E24:S37','Sheet',sheet);
	
	P_weighted = P_matrix .* JPD / sum(JPD,'all');
	P_elec = sum(P_weighted(:));

	% coefficient of variance (normalized standard deviation) of power
	P_var = std(P_matrix(:), JPD(:)) / P_elec;
	c_v = P_var * 100; % convert to percentage

end