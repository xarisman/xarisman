function [result] = g_erotima_iii(p)
    a = p(1); tc = p(2); tw = p(3);
    a1 = a; a2 = a; a3 = a;
    tc_low = tc; tc_upp = tc;
    phi = pi:0.01*pi:5*pi;
    T = zeros(1, length(phi));
    counter1 = 1;
    
    for jj=phi
        ph1 = jj;
        ph2 = jj;
        ph3 = jj;
    
        z_upp = (tc_upp-a3*exp(-(1i)*ph3))/(1-tc_upp*a3*exp(-(1i)*ph3));
        z_low = (tc_low-a2*z_upp*exp(-(1i)*ph2))/(1-tc_low*a2*z_upp*exp(-(1i)*ph2));
        
        T_num = tw^2 + (a1^2)*(abs(z_low)^2)-2*a1*abs(z_low)*tw*cos(ph1-angle(z_low));
        T_denum = 1+(a1^2)*(abs(z_low)^2)*(tw^2)-2*a1*abs(z_low)*tw*cos(ph1-angle(z_low));
            
        T(counter1) = T_num/T_denum;
        
    
        if 10*log(abs(T(counter1)))<-40
            T(counter1) = NaN;
        end
        counter1=counter1+1;
    end
    
    figure;
    hold on;
    plot(phi, (T), 'LineWidth', 2, 'Color', [0,1,0])
    hold off;
    title("a="+ a + "    tc=" + tc + "     tw=" +tw)
end