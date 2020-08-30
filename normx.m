function [y,ic]=normx(x,xmin,xmax,ymin,ymax)
% x is a matrix
% xmin is set as a xmin, but it may not equal to the minimum of x
% xmax the same
% y is the normalized data
% ymin correspond to xmin
% ymax the same
% ic is a notation that x is out of the range of [xmin, xmax]
y = (ymax - ymin)*(x - xmin)/(xmax - xmin) + ymin;
ic = min(x)>=xmin && max(x)<=xmax;


