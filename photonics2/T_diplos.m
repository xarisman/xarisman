function [Td] = T_diplos(twi, tci, a1_i, a2_i, start_phi, end_phi, num_points, df)
    tw = twi;
    tc = tci;
    a1 = a1_i;
    a2 = a2_i;
    phi = linspace(start_phi, end_phi, num_points);
    counter1=1;
    T = zeros(1, length(phi));

    for jj=phi
        ph1 = jj;
        ph2 = jj+df;
        zhta = (tc-a2*exp(-(1i)*ph2))/(1-tc*a2*exp(-(1i)*ph2));
        T_num = tw^2 + (a1^2)*(abs(zhta)^2)-2*a1*abs(zhta)*tw*cos(ph1-angle(zhta));
        T_denum = 1+(a1^2)*(abs(zhta)^2)*(tw^2)-2*a1*abs(zhta)*tw*cos(ph1-angle(zhta));
        
        T(counter1) = T_num/T_denum;
        counter1=counter1+1;
    end
    
    Td = T(:);
end