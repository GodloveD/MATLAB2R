function [ss,r,f] = hyperb(p,x,y)  % hyperb = function to be fitted
vmax = p(1);   km = p(2);          %   by adjusting p = [ vmax, km ].
f = vmax*x./(km+x);                % f = fitting function.
r = y-f;                           % r = residuals.
ss = sum(r.*r);                    % ss = sum of squares.