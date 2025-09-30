clear;close all

p = parameters();
p.control_type = 'reactive';
check_steepness(p)
p.control_type = 'damping';
check_steepness(p)

function check_steepness(p)
    % steepness
    A = p.Hs / (2*sqrt(2));
    w = 2*pi./p.T;
    k = dispersion(w,p.h,p.g);
    steepness = k.*A;
    steepness(p.JPD==0) = NaN;
    
    figure
    mycontour(p,steepness,'Wave steepness $\epsilon$')
    
    % nominal design amplitude to draft ratio
    b = var_bounds();
    X = [b.X_noms;1];
    [~,~,~,val] = simulation(X,p);
    
    T_f_2 = b.T_f_2_nom;
    T_s = p.T_s_over_D_s * b.D_s_nom;
    
    float_amp_draft_ratio = val.X_f / (p.h - T_f_2);
    spar_amp_draft_ratio  = val.X_s / (p.h - T_s);
    
    figure
    subplot 121
    mycontour(p,float_amp_draft_ratio,['Float amplitude-draft ratio $X_f/(h-T_{f,2})$ - ' p.control_type ' control'])
    
    subplot 122
    mycontour(p,spar_amp_draft_ratio,['Spar amplitude-draft ratio $X_s/(h-T_s)$ - ' p.control_type ' control'])
    
    % ratio of amplitude-draft ratio to steepness
    figure
    subplot 121
    mycontour(p, float_amp_draft_ratio ./ steepness,...
        ['Float amplitude indicator $\frac{X_f/(h-T_{f,2})}{\epsilon}$ - ' p.control_type ' control'],[0 1.2]);
    
    subplot 122
    mycontour(p, spar_amp_draft_ratio ./ steepness,...
        ['Spar amplitude indicator $\frac{X_s/(h-T_s)}{\epsilon}$ - ' p.control_type ' control'],[0 1.2]);
end

function mycontour(p,var,mytitle,clims)
    if nargin<4
        levels = 10;
    else
        levels = [min(var(:)), max(var(:)), linspace(clims(1),clims(2),10)];
    end
    contourf(p.T,p.Hs,var,sort(levels))
    xlabel('Wave period T')
    ylabel('Significant wave height H_s')
    title(mytitle,'Interpreter','latex')
    colorbar
    if nargin==4
        clim(clims)
    end
    improvePlot
end