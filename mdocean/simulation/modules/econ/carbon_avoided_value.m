function [avoided_dollars, scc_per_year, avoided_dollars_per_year] = carbon_avoided_value(years, avoided_co2_tons, scc_case)
% carbon_avoided_value
% years: vector of years (e.g., 2025:2050)
% avoided_co2_tons: vector, same length, avoided CO2 each year (metric tons CO2/yr)
% scc_case: "5pct" | "3pct" | "2p5pct" | "high"
%
% Returns:
% avoided_dollars          total avoided $ over all years (SCC table dollars)
% scc_per_year       SCC used each year ($/metric ton CO2)
% avoided_dollars_per_year avoided $ each year

    % SCC table (2007 $ / metric ton CO2)
    scc_years  = [2010 2015 2020 2025 2030 2035 2040 2045 2050];
    scc_5pct   = [10 11 12 14 16 18 21 23 26];
    scc_3pct   = [31 36 42 46 50 55 60 64 69];
    scc_2p5pct = [50 56 62 68 73 78 84 89 95];
    scc_high   = [86 105 123 138 152 168 183 197 212];

    switch string(scc_case)
        case "5pct",   scc_vals = scc_5pct;
        case "3pct",   scc_vals = scc_3pct;
        case "2p5pct", scc_vals = scc_2p5pct;
        case "high",   scc_vals = scc_high;
        otherwise, error("scc_case must be: '5pct','3pct','2p5pct','high'");
    end

    years = years(:);
    avoided_co2_tons = avoided_co2_tons(:);
    if numel(years) ~= numel(avoided_co2_tons)
        error("years and avoided_co2_tons must have the same length.");
    end

    % error if past 2050

    % 1) interpolate SCC for each year
    scc_per_year = interp1(scc_years, scc_vals, years_clamped, "linear");

    % 2) multiply by avoided CO2 each year
    avoided_dollars_per_year = scc_per_year .* avoided_co2_tons;

    % 3) sum across years
    avoided_dollars = sum(avoided_dollars_per_year);
end