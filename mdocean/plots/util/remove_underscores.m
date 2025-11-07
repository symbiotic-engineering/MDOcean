function [cell_str_out] = remove_underscores(cell_str_in)
    % REMOVE_UNDERSCORES Replaces underscores with spaces and capitalizes each word

    if iscell(cell_str_in)
        cell_str_out = cellfun(@remove_underscores_str, cell_str_in, 'UniformOutput', false);
    elseif ischar(cell_str_in)
        cell_str_out = remove_underscores_str(cell_str_in);
    end

end

function [str_out] = remove_underscores_str(str_in)
    % REMOVE_UNDERSCORES Replaces underscores with spaces and capitalizes each word
    cell_split = strsplit(str_in, '_');
    cell_capital = cellfun(@mlreportgen.utils.capitalizeFirstChar, cell_split, 'UniformOutput', false);
    cell_spaces = strcat(cell_capital, " ");
    str_out = horzcat(cell_spaces{:});
end
