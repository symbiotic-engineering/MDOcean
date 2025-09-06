function write_vars_to_tex(varargin)

    % Open the file for writing
    fid = fopen('pubs/journal-paper/numbers.tex', 'w');

    % Write the LaTeX commands for each variable
    for i = 1:nargin
        fprintf(fid, '\\newcommand{\\%s}{%g}\n', inputname(i), varargin{i});
    end

    % Close the file
    fclose(fid);
end