function RM3 = validation_inputs()

RM3 = struct(   'mass_f',       208e3, ...  % p159 RM3 report
                'mass_vc',      224e3, ...  % p159 RM3 report
                'mass_rp',      245e3, ...  % p159 RM3 report
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
                'FOS_b',        3 );        % p158 RM3 report

end