clear all; 
%clc;

N = 50; M=50; dx = 0.2*10^-6; dy = 0.2*10^-6; dz = 0.2*10^-6; lamda = 10^-6;  k_0 = 2*pi/lamda; k_ref = k_0;
n =2; m =2; n_i_j=1;
a = 0.501; w_0 = 1*10^-6;


rx = -a/(dx^2);
ry = -a/(dy^2);
rxx = (1+a)/dx^2;
ryy = (1+a)/dy^2;
C = 2*1j*k_ref/dz+2*a/dx^2+2*a/dy^2-a*(k_0^2*n_i_j^2-k_ref^2);
D = 2*1j*k_ref/dz-2*(1-a)/dx^2+2*(1-a)/dy^2-(1-a)*(k_0^2*n_i_j^2-k_ref^2);

S = C*diag(ones(M*M,1))+ry*diag(ones(M*M-1,1),1)+ry*diag(ones(M*M-1,1),-1)+rx*diag(ones(M*M-N,1),N)+rx*diag(ones(M*M-N,1),-N);
B = D*diag(ones(M*M,1))-ry*diag(ones(M*M-1,1),1)-ry*diag(ones(M*M-1,1),-1)-rx*diag(ones(M*M-N,1),N)-rx*diag(ones(M*M-N,1),-N);
Ss = sparse(S);
Bs = sparse(B);

E_0 = zeros(N-1, M-1);

%E_0(x, y, :) = exp(-(x^2+y^2)/w_0^2);
for ii=0:N-1
    for jj =0:M-1
        E_0(ii+1, jj+1) = exp(-((dx*(ii-25))^2+(dy*(jj-25))^2)/w_0^2);
        %E_0(ii+1, jj+1) = exp(-((ii*10^-6-25*10^-6)^2+(jj*10^-6-25*10^-6)^2)/w_0^2);
        %disp(exp(-((i*10^-6-25*10^-6)^2+(j*10^-6-25*10^-6)^2)/w_0^2))
    end
end

%{
%E_0 = sparse(E_0);
r=1;
rx = a/(dx^2);
ry = a/(dy^2);
rz = (2*i*k_0/dz-2*rx-2*ry);
%rz = 2*1i*k_ref/(dz^2) - a(k_0^2*z-k_ref^2);

% Population of tridiagonal matrix
%A = 2*(1+r)*diag(ones(M-1,1))-r*diag(ones(M-2,1),1)-r*diag(ones(M-2,1),-1);
%A = rx*diag(ones(M-3,1),2)+rx*diag(ones(M-3,1),-2) + ry*diag(ones(M-2,1), 1) + ry*diag(ones(M-2,1), -1) + rz*diag(ones(M-1,-1));


%figure; spy(E_0);
%}

% E_n = reshape(E_0, [], 1);
% E_next = zeros(M^2, 1);
% E_next2 = zeros(M^2, 1); 
% 
% E = zeros(0, 0);
% E_next = Ss\Bs*E_n;
% counter =0;
% rr2 = Ss\Bs;
% 
% for ii=1:M
%     for jj=1:M
%         p = (jj)+counter;
%         E_next2(p) = rr2(p,:)*E_next2;
%     end
%     counter = (jj)+counter;    
% end
% 
% %E2 = abs(E_next-E_next2);
% E = [E, E_next2];



% for i=0:60
%     counter =0;
%     for ii=1:M
%         for jj=1:M
%             p = (jj)+counter;
%             E_next2(p) = rr2(p,:)*E_next2;
%         end
%         counter = (jj)+counter;    
%     end
%     E_next2 = rr2*E_next2;
%     E = [E, E_next2];
% end
E_n = reshape(E_0, [], 1);

E_next = zeros(M^2, 1);  
%disp((E(n, m)+E(n-1,m)-2*E(n, m))/(dx^2))
E = zeros(0, 0);
E_next = Ss\Bs*E_n;
E_prev = E_next;
E = [E, E_next];
%rr = inv(Ss)*Bs;
rr2 = Ss\Bs;
%spy(rr2)
for i=0:60
    E_next = rr2*E_next;
    E = [E, E_next];
end

index = 50;
% Plot the contour
selected_columns = [index];
selected_values = E(:, index);
%selected_values = reshape(selected_values, 1);
%{
Array = zeros(M, N);
counter = 0;
counter2 = 1;
for ii=0:N
    disp(ii)
    while 1==1
        if(counter<51)
            Array(counter+1, ii+1) = selected_values(counter2);
            counter2 = counter2+1
            counter = counter+1;
        else
            counter = counter+1;
            break;
        end
    end
end
%}
selected_values = reshape(selected_values, M, N);

% If you want to transpose the result to put the first 50 values in the first row, and so on
%selected_values = transpose(selected_values);
% Plotting the filled contour of selected values
max_values = max(selected_values(:));
selected_values_normalized = selected_values/max_values;
contourf(abs(selected_values_normalized));


z = index*dz;

z_0 = pi*w_0^2/lamda;
%r = (x.^2+y.^2);
w_z = sqrt(w_0^2*(1+lamda^2*z^2/(pi^2*w_0^4)));
R_z = z*(1+pi^2*w_0^4/(lamda^2*z^2));

%E_r_z = w_0/w_z*exp(-r^2/(w_r^2))*exp(-(1j)*k_0.*r^2/(2.*R_z)).*exp(-(1j)*(k_0*z-tan(z./z_0)));

E_r_z = zeros(M, M);
for ii=1:N
    for jj=1:M
        r = sqrt(((dx*(ii-25))^2+(dy*(jj-25))^2));
        E_r_z(ii, jj) = w_0/w_z*exp(-r^2/(w_z^2))*exp(-(1j)*k_0*r^2/(2*R_z))*exp(-(1j)*(k_0*z-atan(z/z_0)));
    end
end
max_values = max(E_r_z(:));
E_r_z_normalized = E_r_z/max_values;

Errors = abs(E_r_z_normalized - selected_values_normalized);
e = mean(Errors(:))