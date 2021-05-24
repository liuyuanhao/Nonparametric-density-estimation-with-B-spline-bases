function [me_int] = ym(tv,bsn_int,order,knott,index_i,index_j)

%% compute the inner product of the i,jth basis function

w = length(tv)-1;
I = zeros(1,w);
pp = index_i-index_j;
for i = 1:w
    a = tv(i);
    b = tv(i+1);

    
    A = ((a+b)/2-sqrt(3/5)*((b-a)/2));
    B = ((a+b)/2);
    C = ((a+b)/2+sqrt(3/5)*((b-a)/2));
    a_0 = zeros(1,bsn_int);
    a_1 = zeros(1,bsn_int);
    a_2 = zeros(1,bsn_int);
    

    a_0=bspline_basis(index_j-1, order, knott, A)*bspline_basis(index_i-1, order, knott, A);
    a_1=bspline_basis(index_j-1, order, knott, B)*bspline_basis(index_i-1, order, knott, B);
    a_2=bspline_basis(index_j-1, order, knott, C)*bspline_basis(index_i-1, order, knott, C);

    a00=a_0;
    a01=a_1;
    a02=a_2;
    I(1,i) = ((b-a)/18)*(5*a00+8*a01+5*a02);
end
me_int = sum(I);