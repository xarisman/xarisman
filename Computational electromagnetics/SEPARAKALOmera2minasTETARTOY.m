clear all; clc;

N = 50; M=50; dx = 0.2*10^-6; dy = 0.2*10^-6; dz = 0.1*10^-6; lamda = 10^-6;  k_0 = 2*pi/lamda; k_ref = 1; z = 1; n =2; m =2;
a = 0.51; w_0 = 1*10^-6;
%U = zeros(N+1,M+1); r = dt/(a*dx^2);
%U(1,:) = T0*sin((0:M)*pi/M);        % Temperature at t=0 is sinusoidal
%E(N, M) = 0;
E_0 = zeros(N-1, M-1);

%E_0(x, y, :) = exp(-(x^2+y^2)/w_0^2);
for i=0:N-1
    for j =0:M-1
        E_0(i+1, j+1) = exp(-(i^2+j^2)/w_0^2);
    end
end

%E_0 = sparse(E_0);
r=1;
rx = a/(dx^2);
ry = a/(dy^2);
rz = (2*i*k_0/dz-2*rx-2*ry);
%rz = 2*1i*k_ref/(dz^2) - a(k_0^2*z-k_ref^2);

% Population of tridiagonal matrix
%A = 2*(1+r)*diag(ones(M-1,1))-r*diag(ones(M-2,1),1)-r*diag(ones(M-2,1),-1);
%A = rx*diag(ones(M-3,1),2)+rx*diag(ones(M-3,1),-2) + ry*diag(ones(M-2,1), 1) + ry*diag(ones(M-2,1), -1) + rz*diag(ones(M-1,-1));
S = rz*diag(ones(M*M,1))+ry*diag(ones(M*M-1,1),1)+ry*diag(ones(M*M-1,1),-1)+rx*diag(ones(M*M-50,1),50)+rx*diag(ones(M*M-50,1),-50);
B = (1+2*rx+2*ry)*diag(ones(M*M,1))-ry*diag(ones(M*M-1,1),1)-ry*diag(ones(M*M-1,1),-1)-rx*diag(ones(M*M-50,1),50)-rx*diag(ones(M*M-50,1),-50);
Ss = sparse(S);
Bs = sparse(B);

%figure; spy(E_0);
E_n = reshape(E_0, [], 1);

E_next = zeros(M^2, 1);  
%disp((E(n, m)+E(n-1,m)-2*E(n, m))/(dx^2))
E = zeros(0, 0);
E_next = inv(Ss)*Bs*E_n;
E_prev = E_next;
E = [E, E_next];
for i=0:10
    E_next = inv(Ss)*Bs*E_next;
    E = [E, E_next];
end
%{
for n = 2 : N-1           % steps in time
    for m = 2 : M-1       % population of RHS
        %{
        %B(m-1) = r*U(n,m+1)+2*(1-r)*U(n,m)+r*U(n,m-1);
        %B(n-1, m-1) = (E(n, m)+E(n-1, m)-2*E(n, m))/(dx^2);
        %B = rx*diag(ones(M-1,1))+ry*diag(ones(M-2,1),1)+rz*diag(ones(M-2,1),-1)+ry*diag(ones(M-3,1),2)+ry*diag(ones(M-3,1),-2);
        %B(n-1, m-1) = E_0(n-1, m-1);
        %B = (1+2*rx+2*ry)*diag(ones(M*M,1))-ry*diag(ones(M*M-1,1),1)-ry*diag(ones(M*M-1,1),-1)-rx*diag(ones(M*M-4,1),4)-rx*diag(ones(M*M-4,1),-4);
        %}
        

    end
    %X = A\B;             % solution of tridiagonal system
    %E(2:N,2:M) = X';     % storage of solution at new position
end
%}