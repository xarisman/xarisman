clear all; %clc;
tic
Lx = 10*10^-6; Ly = Lx;
N = 30; M=N; dx = Lx/M; dy = Ly/N; dz = 0.04*10^-6; 
lamda = 1.310*10^-6; 
k_0 = 2*pi/lamda; 
R = 5*10^-6; 
a = 0.5; w_0 = 1*10^-6;

n1 = 0; n2 = 1.45; n3 = (n1+n2)/2;
n_i_j = zeros(N, N);
k_ref = k_0*n2;

for ii=1:N
    for jj=1:N
        if sqrt((ii-N/2)^2+(jj-N/2)^2)<R/Lx*N
            n_i_j(ii, jj) = n1;
        elseif sqrt((ii-N/2)^2+(jj-N/2)^2)==R/Lx*N
            n_i_j(ii, jj) = n3;
        else
            n_i_j(ii, jj) = n2;
        end
     end
end


E_0 = zeros(N-1, M-1);
for ii=0:N-1
    for jj =0:M-1
        E_0(ii+1, jj+1) = exp(-((dx*(ii-N/2))^2+(dy*(jj-N/2))^2)/w_0^2);
    end
end

rx = -a/(dx^2);
ry = -a/(dy^2);
rxx = (1-a)/dx^2;
ryy = (1-a)/dy^2;
C = 2*1j*k_ref/dz+2*a/dx^2+2*a/dy^2-a*(k_0^2*n_i_j^2-k_ref^2);
D = 2*1j*k_ref/dz-2*(1-a)/dx^2-2*(1-a)/dy^2+(1-a)*(k_0^2*n_i_j^2-k_ref^2);

S = C.*diag(ones(M*M,1))+ry*diag(ones(M*M-1,1),1)+ry*diag(ones(M*M-1,1),-1)+rx*diag(ones(M*M-N,1),N)+rx*diag(ones(M*M-N,1),-N);
B = D.*diag(ones(M*M,1))+ryy*diag(ones(M*M-1,1),1)+ryy*diag(ones(M*M-1,1),-1)+rxx*diag(ones(M*M-N,1),N)+rxx*diag(ones(M*M-N,1),-N);
Ss = sparse(S);
Bs = sparse(B);
E_n = reshape(E_0, [], 1);

E_next = zeros(M^2, 1);  

E = zeros(0, 0);
E_next = Ss\Bs*E_n;


E = [E, E_next];
rr2 = Ss\Bs;

Gamma = 0;

for ii=0:(12*(10^-4)/dz)
    E_next = rr2*E_next;
    E = [E, E_next];
    if(abs(E_next*exp(k_ref*ii*dz))==inf)
        Gamma = ii*dz;
        disp("aaaaaaaa")
        break
    end
    %disp(ii)
end

%Gamma+dz
E_next = rr2*E_next;
E = [E, E_next];

% clear all; %clc;
% tic
% N = 100; M=N; dx = 10^-5/N; dy = 10^-5/M; dz = 0.04*10^-6; lamda = 10^-6;  k_0 = 2*pi/lamda; k_ref = k_0;
% n_i_j=1;
% a = 0.501; w_0 = 1*10^-6;
E_gamma_dz = E(:, end);
E_gamma = E(:, end-1);

E_gamma_dz = reshape(E_gamma_dz, [], 1);
E_gamma = reshape(E_gamma, [], 1);
%b_01 = kref +1/dz*ln(1);

double_sum = 0;

for ii=0:dx:Lx*N
    for jj=0:dy:Ly*N
        double_sum = double_sum + E_gamma(ii, jj)
    end
end
b_0 = k_ref + log(double_sum)/dz;



