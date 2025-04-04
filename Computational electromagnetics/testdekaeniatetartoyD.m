clear all; %clc;
tic
Lx = 6*10^-6; Ly = Lx;
hx = 10^-6; wx = 10^-6;
N = 111; M=N; dx = Lx/N; dy = Ly/M; dz = 0.2*10^-6; 
lamda = 1.550*10^-6; 
k_0 = 2*pi/lamda; 
%R = 5*10^-6; 
a = 0.50105; w_0 = 1*10^-6;

n1 = 1.6; n2 = 1.45; n3 = 1;
n4 = (n1+n2)/2; n5 = (n1+n3)/2; n6 = (n3+n2)/2;
n7 = (2*n2+n1+n3)/4; n8 = (3*n3+n2)/4;

n_i_j = zeros(N, N);
k_ref = k_0*n4; %ν4

for ii=1:N
    for jj=1:N
        if abs(ii-N/2)<wx/dx && abs(jj-N/2)<hx/dy
            if ii-N/2==-wx/dx+0 || abs(jj-N/2)==hx/dx-0
                if (jj - N/2 == -hx/dy+1 || jj - N/2 == hx/dy-1) && ii-N/2 == -wx/dx+1
                    n_i_j(ii, jj) = n8;
                else
                    n_i_j(ii, jj) =n5;
                end
            else
                n_i_j(ii, jj) = n1;
            end
        elseif (ii-N/2) == wx/dy
            if abs(jj-N/2) < hx/dy-0
                n_i_j(ii, jj) = n4;
            elseif abs(jj-N/2) == hx/dy-0
                n_i_j(ii, jj) = n7;
            else
                 n_i_j(ii, jj) = n6;
            end
        elseif (ii-N/2)> wx/dy+1 %0      
            n_i_j(ii, jj) = n2;
        else
            n_i_j(ii, jj) = n3;
        end
     end
end

tic
E_0 = zeros(N-1, M-1);
for ii=0:N-1
    for jj =0:M-1
        E_0(ii+1, jj+1) = exp(-((dx*(ii-N/2))^2+(dy*(jj-N/2))^2)/w_0^2);
    end
end
toc
rx = -a/(dx^2);
ry = -a/(dy^2);
rxx = (1-a)/dx^2;
ryy = (1-a)/dy^2;
C = 2*1j*k_ref/dz+2*a/dx^2+2*a/dy^2-a*(k_0^2*n_i_j^2-k_ref^2);
D = 2*1j*k_ref/dz-2*(1-a)/dx^2-2*(1-a)/dy^2+(1-a)*(k_0^2*n_i_j^2-k_ref^2);

C = C(:);
D = D(:);

S = zeros(M^2, M^2);
for ii = 1:M^2-1
    for jj = 1:M^2-1
        if(ii<M+1)
            %disp(ii+1)
            S(ii, ii) = C(ii);
            S(ii, ii+1) = ry;
            if ii~=1
                S(ii, ii-1) = ry;
            end
            S(ii, M+ii) = rx;
        elseif(ii<M^2-M+1)
            %disp(ii)
            %disp(ii-M+1)
            S(ii, ii) = C(ii);
            S(ii, ii+1) = ry;
            S(ii, ii-1) = ry;
            S(ii, M+ii) = rx;
            S(ii, -M+ii) = rx;
        else
            S(ii, ii) = C(ii);
            S(ii, ii+1) = ry;
            S(ii, ii-1) = ry;
            S(ii, -M+ii) = rx;
        end
    end
end

B = zeros(M^2, M^2);
for ii = 1:M^2-1
    for jj = 1:M^2-1
        if(ii<M+1)
            %disp(ii+1)
            B(ii, ii) = D(ii);
            B(ii, ii+1) = ryy;
            if ii~=1
                B(ii, ii-1) = ryy;
            end
            B(ii, M+ii) = rxx;
        elseif(ii<M^2-M+1)
            %disp(ii)
            %disp(ii-M+1)
            B(ii, ii) = D(ii);
            B(ii, ii+1) = ryy;
            B(ii, ii-1) = ryy;
            B(ii, M+ii) = rxx;
            B(ii, -M+ii) = rxx;
        else
            B(ii, ii) = D(ii);
            B(ii, ii+1) = ryy;
            B(ii, ii-1) = ryy;
            B(ii, -M+ii) = rxx;
        end
    end
end
toc
S(N^2, N^2) = C(N);
S(N^2, N^2-1) = ry;
S(N^2, -M+N^2) = rx;
B(N^2, N^2) = D(N);
B(N^2, N^2-1) = ryy;
B(N^2, -M+N^2) = rxx;
%Ss = sparse(S);
%spy(Ss)
%Bs = sparse(B);
%spy(Bs)


% S = C.*diag(ones(M*M,1))+ry*diag(ones(M*M-1,1),1)+ry*diag(ones(M*M-1,1),-1)+rx*diag(ones(M*M-N,1),N)+rx*diag(ones(M*M-N,1),-N);
% B = D.*diag(ones(M*M,1))+ryy*diag(ones(M*M-1,1),1)+ryy*diag(ones(M*M-1,1),-1)+rxx*diag(ones(M*M-N,1),N)+rxx*diag(ones(M*M-N,1),-N);

Ss = sparse(S);
Bs = sparse(B);
E_n = reshape(E_0, [], 1); 
E_next = zeros(M^2, 1);  

E = zeros(0, 0);
E_next = Ss\Bs*E_n;
E = gpuArray([E, E_next]); %447secs/431
rr2 = Ss\Bs;
toc
Gamma = 0;
% 
for ii=0:(20*(10^-4)/dz)
    E_next = rr2*E_next;
    E = [E, E_next];
    if(abs(E_next*exp(k_ref*ii*dz))==inf)
        Gamma = (ii-1)*dz;
        disp("aaaaaaaa")
        break
    end
    %disp(abs(E_next*exp(k_ref*ii*dz)))
end

E_gamma_dz = E(:, end);
E_gamma = E(:, end-1);
% 
E_gamma_dz = reshape(E_gamma_dz, [], N);
E_gamma = reshape(E_gamma, [], N);

%B_0_1
double_sum1 = 0;
double_sum2 = 0;
toc 
for ii=1:N
     for jj=1:N
        double_sum1 = double_sum1 + E_gamma(ii, jj);
     end
end
for ii=1:N
     for jj=1:N
        double_sum2 = double_sum2 + E_gamma_dz(ii, jj);
     end
end 

b_0_1 = k_ref + log(double_sum2/double_sum1)/dz
toc
%B_0_2
% %B_0_2
 double_sum_12 = 0;
 double_sum_22 = 0;

 for ii=1:N
     for jj=1:N
        double_sum_12 = double_sum_12 + log(abs(E_gamma_dz(ii, jj))/abs(E_gamma(ii,jj)))*(abs(E_gamma(ii,jj))^2);
     end
end
for ii=1:N
     for jj=1:N
        double_sum_22 = double_sum_22 + abs(E_gamma(ii, jj))^2;
     end
end

b_0_2 = k_ref + 1/dz*double_sum_12/double_sum_22

%b_0_3
U_0 = E_gamma_dz;
double_sum_13 = 0;
double_sum_23 = 0;

for ii=1:N
     for jj=1:N
         %AA = diff(E_0, 2);
         %A = abs(diff(diff(U_0, N-2), N));
         A = abs(diff(U_0, 1));
        double_sum_13 = double_sum_13  + (k_0^2*n_i_j(ii, jj)^2-k_ref^2)*abs(U_0(ii, jj))^2-A(ii)^2;
     end
end
for ii=1:N
     for jj=1:N
        double_sum_23 = double_sum_23 + 2*k_ref*abs(U_0(ii, jj))^2;
     end
end

b_0_3 = k_ref + double_sum_13/double_sum_23
N1 = b_0_1/k_0;
N2 = b_0_2/k_0;
N3 = b_0_3/k_0;
toc