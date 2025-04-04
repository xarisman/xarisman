clear all; %clc;
tic
N = 100; M=N; dx = 10^-5/N; dy = 10^-5/M; dz = 0.04*10^-6; lamda = 10^-6;  k_0 = 2*pi/lamda; k_ref = k_0;
n_i_j=1;
a = 0.501; w_0 = 1*10^-6;


rx = -a/(dx^2);
ry = -a/(dy^2);
rxx = (1-a)/dx^2;
ryy = (1-a)/dy^2;
C = 2*1j*k_ref/dz+2*a/dx^2+2*a/dy^2-a*(k_0^2*n_i_j^2-k_ref^2);
D = 2*1j*k_ref/dz-2*(1-a)/dx^2-2*(1-a)/dy^2+(1-a)*(k_0^2*n_i_j^2-k_ref^2);

S = C*diag(ones(M*M,1))+ry*diag(ones(M*M-1,1),1)+ry*diag(ones(M*M-1,1),-1)+rx*diag(ones(M*M-N,1),N)+rx*diag(ones(M*M-N,1),-N);
B = D*diag(ones(M*M,1))+ryy*diag(ones(M*M-1,1),1)+ryy*diag(ones(M*M-1,1),-1)+rxx*diag(ones(M*M-N,1),N)+rxx*diag(ones(M*M-N,1),-N);
Ss = sparse(S);
Bs = sparse(B);

E_0 = zeros(N-1, M-1);
for ii=0:N-1
    for jj =0:M-1
        E_0(ii+1, jj+1) = exp(-((dx*(ii-N/2))^2+(dy*(jj-N/2))^2)/w_0^2);
    end
end
toc
E_n = reshape(E_0, [], 1);

E_next = zeros(M^2, 1);  

E = zeros(0, 0);
E_next = Ss\Bs*E_n;


E = [E, E_next];
rr2 = Ss\Bs;



for i=0:(12*(10^-6)/dz)
    E_next = rr2*E_next;
    E = [E, E_next];
end
toc
%z = 10*10^-6;
for iii = [1, 2, 5, 10, 12]
    z = iii*10^-6;
    index = z/dz;
    
    selected_values = E(:, fix(index));
    selected_values = reshape(selected_values, M, N);
    max_values = max(selected_values(:));
    
    selected_values_normalized = selected_values/max_values;
    figure()
    contourf(abs(selected_values_normalized));
    
    
    
    z_0 = pi*w_0^2/lamda;
    
    w_z = sqrt(w_0^2*(1+lamda^2*z^2/(pi^2*w_0^4)));
    R_z = z*(1+pi^2*w_0^4/(lamda^2*z^2));
    
    E_r_z = zeros(M, M);
    
    for ii=1:N
        for jj=1:M
            r = sqrt(((dx*(ii-N/2))^2+(dy*(jj-N/2))^2));
            E_r_z(ii, jj) = (w_0/w_z)*exp(-r^2/(w_z^2))*exp(-(1j)*k_0*r^2/(2*R_z))*exp(-(1j)*(k_0*z-atan(z/z_0)));
        end
    end
    % max_values = max(E_r_z(:));
    % E_r_z_normalized = E_r_z/max_values;
    
    %Errors = abs(E_r_z_normalized - selected_values_normalized);
    Errors = abs(E_r_z - selected_values);
    e = mean(Errors(:))*M/sum(sum(abs(selected_values)))
end
toc