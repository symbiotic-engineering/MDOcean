function [varargout] = scale_bessels_symbolic(varargin)
    varargout = cell(1,numel(varargin));
    for idx = 1:numel(varargin)
        in = varargin{idx};
        tmp = scale_bessel(in, 'i');
        out = scale_bessel(tmp,'j');
        varargout{idx} = out;
    end
end

function expr_scaled = scale_bessel(expr,ij)
    syms mu Z
    if strcmp(ij,'i')
        scaled_bessel(mu,Z) = besseli(mu,Z)*exp(abs(real(Z)));
        unscaled_bessel(mu,Z) = besseli(mu,Z);
    elseif strcmp(ij,'j')
        scaled_bessel(mu,Z) = besselj(mu,Z)*exp(abs(imag(Z)));
        unscaled_bessel(mu,Z) = besselj(mu,Z);
    else
        error('ij wrong')
    end
    
    % here the function call to besseli/j is intended to be replaced with
    % besseli/j(mu,Z,1) after matlabFunction generates the file
    
    args = get_besselij_args(expr,ij);
    
    % Substitute all calls to besseli with myFunction
    expr_scaled = expr; % Start with the original expression
    for i = 1:size(args,1)
        old = unscaled_bessel(args(i,1), args(i,2));
        new = scaled_bessel(  args(i,1), args(i,2));
        expr_scaled = subs(expr_scaled, old, new);
    end
end

function args = get_besselij_args(expr,ij)
    if strcmp(ij,'i')
        fcn_string = 'besseli';
        fcn = @(mu,Z)besseli(mu,Z);
    elseif strcmp(ij,'j')
        fcn_string = 'besselj';
        fcn = @(mu,Z)besselj(mu,Z);
    else
        error('ij wrong')
    end

    args = sym([]);
    for array_idx = 1:numel(expr)
        parent = expr(array_idx);

        % if the parent doesn't have a bessel, there are no bessel args to find
        if has(parent,fcn_string) 
            
            % check which children have a bessel
            ch = children(parent);
            ch_has = false(length(ch));
            for idx = 1:length(ch)
                child = ch{idx};
                ch_has(idx) = has(child,fcn_string);
            end 
        
            % if the parent has a bessel and none of the children do, the parent probably is
            % the plain bessel. This logic breaks if one of the arguments of the
            % bessel is a bessel itself, but MEEM doesn't have terms like that so
            % I'll ignore that edge case.
            if all(~ch_has)
                if length(ch)==2 && parent==fcn(ch{:})
                    args = [args; ch{:}];
                else 
                    error('confusion: no children have bessel, but parent isnt plain bessel')
                end
            else
            % if the parent has a bessel and so does 1+ of the children, keep
            % expanding those children.
                ch_with_bessel = ch(ch_has);
                for idx=1:length(ch_with_bessel)
                    arg_idx = get_besselij_args(ch_with_bessel{idx}, ij);
                    args = [args; arg_idx];
                end
            end
        end
    end

end