function [minval maxval] = imrange(inpict)
%   [min max] = IMRANGE(INPICT)
%       returns the global min and max pixel values in INPICT
%       
%   INPICT can be a vector or array of any dimension or class

x = inpict(:);
minval = double(min(x));
maxval = double(max(x));

if nargout < 2
  minval = cat(2,minval,maxval);
end

end



