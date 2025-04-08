function inputCell = modify_suffix(inputCell)
    
    for i=1:length(inputCell)
        inputString = inputCell{i};

        if endsWith(inputString, '_s')
            modifiedString = replace(inputString, '_s', '_spar');
        elseif endsWith(inputString, '_f')
            modifiedString = replace(inputString, '_f', '_float');
        elseif endsWith(inputString, '_vc')
            modifiedString = replace(inputString, '_vc', '_vertical_column');
        elseif endsWith(inputString, '_rp')
            modifiedString = replace(inputString, '_rp', '_reaction_plate');
        elseif endsWith(inputString, '_avg')
            modifiedString = replace(inputString, '_avg', '_average');
        elseif endsWith(inputString, '_tot')
            modifiedString = replace(inputString, '_tot', '_total');
        else
            modifiedString = inputString; % No change if it doesn't match
        end

        inputCell{i} = modifiedString;
    end
end