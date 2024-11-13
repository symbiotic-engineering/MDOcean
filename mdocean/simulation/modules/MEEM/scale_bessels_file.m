function scale_bessels_file(filename,folder)

    new_filename = strcat(filename(1:end-3), 'ed'); % change scaling to scaled
    
    % Open the file for reading
    fid = fopen([folder filename '.m'], 'r');
    if fid == -1
        error('Cannot open the file.');
    end
    
    % Read the content of the file
    file_content = fread(fid, '*char')';
    fclose(fid);
    
    % Replace occurrences of besseli(nu,Z) with besseli(nu,Z,1)
    modified_content = regexprep(file_content,     'besseli\((\w+),\s*(\w+)\)', 'besseli($1,$2,1)');

    % Replace occurrences of besselj(nu,Z) with besselj(nu,Z,1)
    modified_content = regexprep(modified_content, 'besselj\((\w+),\s*(\w+)\)', 'besselj($1,$2,1)');
    
    % Replace function header
    modified_content = regexprep(modified_content, filename, new_filename);

    % Open a new file for writing
    fid = fopen([new_filename '.m'], 'w');
    if fid == -1
        error('Cannot create the new file.');
    end
    
    % Write the modified content to the new file
    fwrite(fid, modified_content);
    fclose(fid);

end