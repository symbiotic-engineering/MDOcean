% ----------------------------------------------------------------------- %
% Function table2latex(T, filename) converts a given MATLAB(R) table into %
% a plain .tex file with LaTeX formatting.                                %
%                                                                         %
%   Input parameters:                                                     %
%       - T:        MATLAB(R) table. The table should contain only the    %
%                   following data types: numeric, boolean, char or string.
%                   Avoid including structs or cells.                     %
%       - filename: (Optional) Output path, including the name of the file.
%                   If not specified, the table will be stored in a       %
%                   './table.tex' file.                                   %  
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       LastName = {'Sanchez';'Johnson';'Li';'Diaz';'Brown'};             %
%       Age = [38;43;38;40;49];                                           %
%       Smoker = logical([1;0;1;0;1]);                                    %
%       Height = [71;69;64;67;64];                                        %
%       Weight = [176;163;131;133;119];                                   %
%       T = table(Age,Smoker,Height,Weight);                              %
%       T.Properties.RowNames = LastName;                                 %
%       table2latex(T);                                                   %                                       
% ----------------------------------------------------------------------- %
%   Version: 1.1                                                          %
%   Author:  Victor Martinez Cagigal                                      %
%   Date:    09/10/2018                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
function table2latex(tab, filename, special_col_spec, special_first_row)
    
    % Error detection and default parameters
    if nargin < 2
        filename = 'table.tex';
        fprintf('Output path is not defined. The table will be written in %s.\n', filename); 
    elseif ~ischar(filename)
        error('The output file name must be a string.');
    else
        if ~strcmp(filename(end-3:end), '.tex')
            filename = [filename '.tex'];
        end
    end
    if nargin < 1, error('Not enough parameters.'); end
    if ~istable(tab), error('Input must be a table.'); end
    if nargin<3
        special_col_spec = [];
    end
    if nargin<4
        special_first_row = [];
    end
    
    % Parameters
    n_col = size(tab,2);
    col_spec = [];
    for c = 1:n_col
        col_spec = [col_spec 'l'];
    end
    col_names = strjoin(tab.Properties.VariableNames, ' & ');
    row_names = tab.Properties.RowNames;
    if ~isempty(row_names)
        col_spec = ['l' col_spec]; 
        col_names = ['& ' col_names];
    end
    
    if ~isempty(special_col_spec)
        col_spec = strrep(special_col_spec, '\', '\\');
    end

    % Writing header
    
    table_string = '';
    table_string = append(table_string, ['\\begin{tabular}{' col_spec '}\n']);

    if ~isempty(special_first_row)
        table_string = append(table_string, [strrep(special_first_row, '\', '\\') ' \n']);
    end

    table_string = append(table_string, [col_names '\\\\ \n']);
    table_string = append(table_string, '\\hline \n');


    % Writing the data
    try
        for row = 1:size(tab,1)
            row_data = cell(1,n_col);
            row_data{1,n_col} = [];
            for col = 1:n_col
                value = tab{row,col};
                
                use_percent = contains(tab.Properties.VariableNames{col},'error','IgnoreCase',true);
                value = format_value(value, use_percent);
                
                row_data{1,col} = char(value);
            end
            if ~isempty(row_names)
                row_data = [row_names{row}, row_data];
            end
            row_string = append(strjoin(row_data, ' & '), ' \\\\ \n'); % need 4 \ in row_string and 2 \ in file
            table_string = append(table_string, row_string);
            clear row_data;
        end
    catch err
        msg = ['Cannot generate latex table string. Make sure that table only ' ...
            'contains chars, strings or numeric values.'];
        base_err = MException('MDOcean:table2latex',msg);
        base_err = addCause(base_err,err);
        throw(base_err)
    end
    
    table_string = append(table_string,'\\end{tabular}');

    fileID = fopen(filename, 'w');
    fprintf(fileID, table_string);

    % Closing the file
    fclose(fileID);
end

function value = format_value(value, use_percent)
    if isstruct(value)
        error('Table must not contain structs.');
    end
    while iscell(value)
        value = value{1,1};
    end
    if ismissing(value)
        value = '-';
    end
    if ~isempty(value) && (ischar(value) || isstring(value))
        value = char(value);
        if isstrprop(value(1), 'digit')
            value = str2double(value);
        end
    end
    
    % Format the output
    if isnumeric(value)
        if isinf(value)
            value = '$\infty$';
        elseif use_percent
            value = sprintf('$%.1f\\\\%%%% $',value*100);
        elseif isscalar(value) % use engineering notation
            value = engr_notation(value);
        else % array
            numbers = strjoin(arrayfun(@(x) engr_notation(x),value,'UniformOutput',false), ', ');
            value = ['[ ' numbers ' ]'];
        end
    else
        value = formattedDisplayText(value,"UseTrueFalseForLogical",true);
    end
end

function formatted = engr_notation(value)
    exponent = floor(log10(abs(value))/3) * 3; % Round down to nearest multiple of 3
    mantissa = value / 10.^exponent; % Calculate mantissa
    if exponent==0
        formatted = sprintf('$%.3g $',mantissa);
    else
        formatted = sprintf('$%.3g \\\\cdot 10^{%d}$', mantissa, exponent);
    end
end