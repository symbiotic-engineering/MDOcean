function u = unicode(codepoint)

% codepoint is a character array representing the part after U+ in unicode.
% u is a character that renders the unicode, for i.e. graphics.
% This function handles the required difference between the way matlab handles
% 4-digit codepoints (on the "basic multilingual plane") and 
% 5- and 6-digit codepoints (on the "supplementary planes").
% Example codepoints: '1f7c0' for U+1F7C0 🟀, '1d54f' for U+1D54F 𝕏, '2102' for U+2102 ℂ

    len = length(codepoint);
    dec = hex2dec(codepoint);
    if len==4
        u = char(dec);
    elseif len==5 || len==6
        u = java.lang.Character.toChars(dec).';
    else
        error('Unicode codepoint must be 4-6 digits long')
    end

end