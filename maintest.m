clc; clear;
n_s = 100000;%amount of samples
order = 3;%the order of the basis functions

%% generate samples 
%rng default;
% mu0 = 0;
% sigma0 = 1;
% f = @(x) (1/sqrt(2*pi )*sigma0)*exp(-(x-mu0).^2/(2*sigma0^2))*(x>=-3 & x<=3);

% lambada = 1;
% f = @(x) lambada*exp(-lambada*x)*(x>=0 & x<=3);

% n = 3;
% f = @(x) 0.8*((x.^(n/2-1)*exp(-x/2))/(gamma(n/2)*2^(n/2)))*(x>=0 & x<=8);

% mu0 = 0; 
% sigma0 = 1;
% mu1 = 4; 
% sigma1 = 1;  
% f = @(x) 0.5*(1/sqrt(2*pi )*sigma0)*exp(-(x-mu0).^2/(2*sigma0^2))*(x>=-3 & x<=7)+...
%     0.5*(1/sqrt(2*pi)*sigma1)*exp(-(x-mu1).^2/(2*sigma1^2))*(x>=-3 & x<=7);

 mu0 = 0;
 sigma0 = 1;
 f = @(x) (1/sqrt(2*pi )*sigma0)*exp(-(x-mu0).^2/(2*sigma0^2))*(x>-5 & x<=5);


% sigma0 = 0.2;
% mu0 = 4; 
% sigma1 = 0.9; 
% mu1 = 1; 
% sigma2 = 1.4;
% mu2 = 8;
% f = @(x) 0.6*(1/sqrt(2*pi)*sigma0)*exp(-(x-mu0).^2/(2*sigma0^2))+...
%     0.2*(1/sqrt(2*pi)*sigma1)*exp(-(x-mu1).^2/(2*sigma1^2))+...
%     0.2*(1/sqrt(2*pi)*sigma2)*exp(-(x-mu2).^2/(2*sigma2^2));

%area = integral(f,-5,5);
N = n_s;
x_sample = slicesample(1,N,'pdf',f,'thin',5,'burnin',n_s);%generate sample
x_sample = sort(x_sample);


% pp=0.5;
% Xn = max(x_sample);X1 = min(x_sample);
% l = X1 - pp*(Xn - X1);
% u = Xn + pp*(Xn - X1);
% 
% h = 0.1;%bandwidth
% a = 1/h;

%% generate B spline basis function
[knott,tv,c,d] = nodesequence(x_sample,n_s);
bsn_int = length(knott) - order;%number of basis function


%% plot reconstruct statistical model
n_t = n_s;
% x_t = linspace(0,1,n_t)*(d-c)+c; %use any point lead to error.
%MUST USE initial x_sample
x_t = x_sample;
y_t = zeros(bsn_int, n_t);
for i = 1:bsn_int
    y_t(i,:) = bspline_basis(i-1, order, knott, x_t);%store basis functions 
end


%% Assemble the system 
A = zeros(bsn_int,bsn_int);
theta = zeros(bsn_int,1);
for i = 1:bsn_int
    for j = 1:bsn_int
        if abs(i-j) > order
            A(i,j) = 0;
        else
            A(i,j) = ym(tv,bsn_int,order,knott,i,j);%%Inner product
        end
    end
    theta(i) = sum(y_t(i,:))/n_t;
end

alpha = A\theta;%% compute

%% calculate the value of true p.d.f. f
F_t=zeros(1,n_s);
for i=1:n_s
    x=x_sample(i,1);
    F_t(1,i)=f(x);
end

%% histogram
h0 = max(x_sample)-min(x_sample);
d0 = h0/100;

figure(1) 
[x_exact1,y_exact1] = hist(x_sample,100);
f0 = x_exact1/n_s;
f1 = f0/d0;
hold on;
h = bar(y_exact1,f1,'hist');
h.FaceColor = [.8 .8 1];

%% solution
ff = zeros(bsn_int,n_t);
for i = 1:bsn_int
    ff(i,:) = alpha(i)*y_t(i,:);
end

finalf = zeros(1,n_t);
for i = 1:n_t
    for j = 1:bsn_int
        finalf(i) = finalf(i) + ff(j,i);
    end
end

plot(x_t,finalf);
hold on


%% exact
x_exact = zeros(n_t,1);
y_exact = zeros(n_t,1);
i = 0;
for x = -5:0.001:5
    i = i+1;
    x_exact(i) = x;
    y_exact(i) = f(x); 
end

plot(x_exact,y_exact);
hold on

%% Calculate MSE/MAE/RMAE/R square
Sum1 = 0;
Sum2 = 0;
Sum3 = 0;
Sum4 = 0;

ave = mean(F_t);

for i = 1:n_s
    sub1 = (finalf(1,i)-F_t(1,i))^2;
    sub2 = abs(finalf(1,i)-F_t(1,i));
    sub3 = sub1;
    sub4 = (ave - F_t(1,i))^2;
    Sum1 = Sum1+sub1;
    Sum2 = Sum2+sub2;
    Sum3 = Sum3+sub3;
    Sum4 = Sum4+sub4;
end

MSE_int = Sum1/n_s
RMSE_int = sqrt(MSE_int)
MAE_int = Sum2/n_s
R_square = 1-(Sum3/Sum4)