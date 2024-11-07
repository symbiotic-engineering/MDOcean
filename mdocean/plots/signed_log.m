function signed_log(Z,C,levels,X,Y)
% modified from the following
% https://www.mathworks.com/matlabcentral/answers/1700655-symmetric-diverging-log-color-scale#answer_1380451
    
    % assign default arguments if not provided
    if nargin<2
        C = [];
    end
    if nargin<3
        levels = [];
    end
    if nargin<4
        use_XY = false;
    elseif nargin==5
        use_XY = true;
    else
        error('Wrong number of input arguments')
    end
    
    if isempty(C)
        % controls smallest order of magnitude near zero
        C = log10(min(nonzeros(abs(Z)),[],'all')); 
    end

    % capital Z is the data, lowercase z is the signed-log-transformed data
    z = signedlog10cont(Z,C); 

    % get the maximum value of abs(Z)
    Z_range = max(abs(imrange(Z)));

    % get color range. this needs to be symmetric unless you also make a 
    % custom asymmetric divergent colormap
    % the bluewhitered colormap is smart and adjusts to all pos/neg so need
    % to account for this in crange.
    if all(z>=0,'all')
        c_ends = [0 1];
    elseif all(z<=0,'all')
        c_ends = [-1 0];
    else
        c_ends = [-1 1];
    end
    crange_Z = Z_range * c_ends;
    crange_z = signedlog10cont(crange_Z,C);

    round = false; % control whether to round up colorbar to integer exponent
    if round
        crange_z = c_ends + fix(crange_z);
    end

    % create ticks, ticklabels
    % choose nticks as desired
    if ~isempty(levels)
        % add the maximum abs value to the specified levels to avoid unfilled 
        % color between [-Inf levels(1)]
        pos_tick_Z_values = sort([Z_range levels]);
    else % find levels automatically if not provided
        max_nticks = 10;
        nticks_temp = max(crange_z)-min(crange_z)+1;
        nticks = min(nticks_temp,max_nticks);
        pos_tick_Z_values = logspace(min(crange_z),max(crange_z),nticks);
    end

    % set Ticks and TickLabels
    tick_Z_values = [-flip(pos_tick_Z_values) 0 pos_tick_Z_values]; % make symmetric
    tick_z_values = signedlog10cont(tick_Z_values,C);
    tick_labels = tick_Z_values;

    % plot z - use contourf if X and Y are provided, imagesc if not  
    if use_XY
        contourf(X,Y,z,'LevelList',tick_z_values);
    else
        imagesc(z)
    end

    % add colorbar with tick labels
    cb = colorbar('Ticks',tick_z_values);
    cb.TickLabels = arrayfun(@(x) sprintf('%.1e', x), tick_labels, 'UniformOutput', false);
    % set limits for the caxis 
    caxis(crange_z); % represented in signed log10

    % use a symmetric colormap
    colormap(bluewhitered())

    improvePlot
    axis square
end

function out = signedlog10cont(in,C)
    % modified continuous signed-log
    % see Measurement Science and Technology (Webber, 2012)  
    % allows for negative exponents up to C
    out = sign(in).*(log10(1+abs(in)/(10^C)));
end

