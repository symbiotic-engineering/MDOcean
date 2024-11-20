function [RM3_actual] = validation_inputs(mode)

RM3_report = struct(...
                'mass_f',       208e3, ...  % p159 RM3 report
                'mass_vc',      224e3, ...  % p159 RM3 report
                'mass_rp',      245e3, ...  % p159 RM3 report
                'mass_tot',     680e3, ...  % p159 RM3 report
                'vol_f',        742.044,... % calculated (MDOcean F24 debugging slide 23)
                'vol_s',        1008, ...   % calculated (MDOcean F24 debugging slide 27)
                'capex',        [17e6 61e6, 207e6, 390e6],  ... 
                ...                         % B15:E15 in CBS Report Tables
                'opex',         [1.2e6, 3.3e6, 6.6e6, 9.4e6],... 
                ...                         % B58:E58 in CBS Report Tables
                'LCOE',         [4.48, 1.45, 0.85, 0.76], ... 
                ...                         % B98:F98 in CBS Report Tables
                ...                         % 4.48 is estimated from Fig 5-33 p175 RM3 report
                'power_avg',    85.9e3, ... % S14 in CBS Performance & Economics
                'power_max',    286e3, ...  % S15 in CBS Performance & Economics
                'force_heave',  8500e3,...  % p156 RM3 report
                'FOS_b',        3,...       % p158 RM3 report
				'c_v',			nominal_c_v()); % calculated from data in CBS Performance & Economics

RM3_wecsim = struct(...
                'vol_f',        725.833, ...% submerged volume from WAMIT RM3.out
                'vol_s',        726.755, ...% submerged volume from WAMIT RM3.out
                'CB_f',         1.2929, ... % center of buoyancy below waterline (MDOcean F24 debugging slide 25)
                'CG_f',         0.2833 ...  % center of gravity below waterline (MDOcean F24 debugging slide 25)
              );

if strcmpi(mode,'wecsim')
    RM3_actual = RM3_wecsim;
elseif strcmpi(mode,'report')
    RM3_actual = RM3_report;
else
    error("Invalid validation mode: use 'wecsim' or 'report.'")
end

end