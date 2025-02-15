function [cell_out] = remove_underscores(cell_in)
%REMOVE_UNDERSCORES Replaces underscores with spaces and capitalizes each word

    cell_out = cellfun(@remove_underscores_str,cell_in,'UniformOutput',false);

end

function [str_out] = remove_underscores_str(str_in)
%REMOVE_UNDERSCORES Replaces underscores with spaces and capitalizes each word
    cell_split = strsplit(str_in,'_');
    cell_capital = cellfun(@mlreportgen.utils.capitalizeFirstChar,cell_split,'UniformOutput',false);
    cell_spaces = strcat(cell_capital," ");
    str_out = horzcat(cell_spaces{:});
end

