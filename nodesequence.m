function [knott,tv,c,d]=nodesequence(x_sample,n_s)
%Input arguments:
% x_sample:
%    random sample
% n_s:
%   the number of sample
%
% Output arguments:
% knott:
%   the initial node sequence
% tv:
%   node vector before extended nodes are added
% c:
%   the interval minimum
% d:
%   the interval maximum

x_sample=sort(x_sample);
c=x_sample(1);
d=x_sample(n_s);
tv = linspace(0,1,50)*(d-c)+c;
knott=[c c tv d d];