%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.


function [x] = minmaxn(x)
a = max(x);
b = min(x);
x = (x - b) ./ (a - b);
end