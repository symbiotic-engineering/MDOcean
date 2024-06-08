function signed_log(A)
% https://www.mathworks.com/matlabcentral/answers/1700655-symmetric-diverging-log-color-scale#answer_1380451
    figure
    % plot A
    % i'm going to plot xdata just for sake of giving some scale information
    imagesc(signedlog10(A)) % represented in signed log10
    % use a symmetric colormap
    colormap(bluewhitered())
    % get the minimum and maximum value of A
    % this needs to be symmetric unless you also make a 
    % custom asymmetric divergent colormap
    % the bluewhitered colormap is smart and adjusts to all pos/neg so need
    % to account for this in crange.

    arange = max(abs(imrange(A)));
    if all(A>=0,'all')
        crange = [0 arange];
        crange = [0 1] + fix(signedlog10(crange));
    elseif all(A<=0,'all')
        crange = [-arange 0];
        crange = [-1 0] + fix(signedlog10(crange));
    else
        crange = [-1 1]*arange; % make symmetric
        crange = [-1 1] + fix(signedlog10(crange));
    end
     % you may choose to round to integer values
    % set limits for the caxis 
    caxis(crange); % represented in signed log10
    % create ticks, ticklabels
    % choose nticks as desired 
    nticks = max(crange)+1;
    ticklabels = logspace(0,max(crange),nticks);
    ticklabels = [-flip(ticklabels) 0 ticklabels]; % make symmetric
    % set Ticks and TickLabels
    colorbar('Ticks',signedlog10(ticklabels),'TickLabels',ticklabels)
    improvePlot
    axis square
end

function out = signedlog10(in)
    % naive signed-log
    %out = sign(in).*log10(abs(in));
    
    % modified continuous signed-log
    % see Measurement Science and Technology (Webber, 2012)
    C = 0; % controls smallest order of magnitude near zero
    out = sign(in).*(log10(1+abs(in)/(10^C)));
end
