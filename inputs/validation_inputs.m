function RM3 = validation_inputs()

RM3 = struct(   'mass_f',       208e3, ...
                'mass_vc',      224e3, ...
                'mass_rp',      245e3, ...
                'capex',        17e6,  ...
                'opex',         1.2e6, ...
                'LCOE',         0.75,  ...
                'power_avg',    286e3 * 0.3, ...
                'power_max',    286e3, ...
                'force_heave',  8500e3,...
                'FOS_b',        3 );

end