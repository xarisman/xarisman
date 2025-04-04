Z_0 = 120*pi; N=2;
n_h = 2.14; h_h = Z_0/n_h;
n_l = 1.46; h_l = Z_0/n_l;
n_g = 1.52; h_g = Z_0/n_g;
n_a = 1; h_a = Z_0/n_a;
counter = 1;
lamda_0 = 860*10^-9;

lambda_values = 400*(10^-9):1*(10^-9):1100*(10^-9);
Reflective_coefficient = zeros(length(lambda_values),1);
T = zeros(length(lambda_values),1);



for lamda = lambda_values
    Z = zeros(N*2+3+1,1);
    h = zeros(N*2+3, 1);
    d = zeros(N*2+3, 1);
    k = zeros(N*2+3, 1);
    
    
    for i=1:(N*2+3)
        if i==1 || i==N*2+3
            d(i) = lamda_0/8/n_l;
        else
            if mod(i, 2)==1
                d(i) = lamda_0/4/n_l;
            else
                d(i) = lamda_0/4/n_h;
            end
        end
    
        if mod(i, 2)==1
            h(i) = h_l;
            k(i) = 2*pi/lamda*n_l;
        else
            h(i) = h_h;
            k(i) = 2*pi/lamda*n_h;
        end
    end
    
    Z(end) = h_g;
    for i=1:(N*2+3)
        Z(end-i) = (Z(end-i+1) + h(end-i+1)*(1i)*tan(k(end-i+1)*d(end-i+1)))/(...
            h(end-i+1)+Z(end-i+1)*(1i)*tan(k(end-i+1)*d(end-i+1)))*h(end-i+1);

    end
    
    Reflective_coefficient(counter) = abs(((Z(1) - h_a)/(Z(1) + h_a)));
    T(counter) = 1-Reflective_coefficient(counter);

    counter = counter+1;
end

plot(lambda_values * 10^9, T);
xlabel('lamda');
ylabel('T');
title('T vs lamda');
grid on;