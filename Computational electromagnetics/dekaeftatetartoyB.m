clear all; %clc;
tic
Lx = 10*10^-6; Ly = Lx;
hx = 10^-6; wx = 10^-6;
N = 100; M=N; dx = Lx/N; dy = Ly/M; dz = 0.4*10^-6; 
lamda = 1.550*10^-6; 
k_0 = 2*pi/lamda; 
R = 5*10^-6; 
a = 0.5; w_0 = 1*10^-6;

n1 = 1.6; n2 = 1.45; n3 = 1;
n4 = (n1+n2)/2; n5 = (n1+n3)/2; n6 = (n3+n2)/2;
n7 = (2*n2+n1+n3)/4; n8 = (3*n3+n2)/4;

n_i_j = zeros(N, N);
k_ref = k_0*n2;

for ii=1:N
    for jj=1:N
        %disp(10^-6/N)
        if abs(ii-N/2)<10^-6/dx && abs(jj-N/2)<10^-6/dy
            if ii-N/2==-10^-6/dx+1 || abs(jj-N/2)==10^-6/dx-1
                if (jj - N/2 == -10^-6/dy+1 || jj - N/2 == 10^-6/dy-1) && ii-N/2 == -10^-6/dx+1
                    n_i_j(ii, jj) = n8;
                else
                    n_i_j(ii, jj) =n5;
                end
            else
                n_i_j(ii, jj) = n1;
            end
        elseif (ii-N/2) == 10^-6/dy
            if abs(jj-N/2) < 10^-6/dy-1
                n_i_j(ii, jj) = n4;
            elseif abs(jj-N/2) == 10^-6/dy-1
                n_i_j(ii, jj) = n7;
            else
                 n_i_j(ii, jj) = n6;
            end
        elseif (ii-N/2)> 10^-6/dy      
            n_i_j(ii, jj) = n2;
        else
            n_i_j(ii, jj) = n3;
        end
     end
end
