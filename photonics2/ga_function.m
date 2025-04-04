function [func] = ga_function(p, BWi, c1i, c2i)
    a = p(1); tc = p(2); tw = p(3);
    a1 = a; a2 = a; a3 = a;
    tc_low = tc; tc_upp = tc;
    phi = pi:0.01*pi:5*pi;
    T = zeros(1, length(phi));
    counter1 = 1;
    counter2 = 1;
    BW = round(BWi/(phi(2)-phi(1)));
    c1 = c1i;
    c2 = c2i;


    for jj=phi
        ph1 = jj;
        ph2 = jj;
        ph3 = jj;
    
        z_upp = (tc_upp-a3*exp(-(1i)*ph3))/(1-tc_upp*a3*exp(-(1i)*ph3));
        z_low = (tc_low-a2*z_upp*exp(-(1i)*ph2))/(1-tc_low*a2*z_upp*exp(-(1i)*ph2));
        
        T_num = tw^2 + (a1^2)*(abs(z_low)^2)-2*a1*abs(z_low)*tw*cos(ph1-angle(z_low));
        T_denum = 1+(a1^2)*(abs(z_low)^2)*(tw^2)-2*a1*abs(z_low)*tw*cos(ph1-angle(z_low));
            
        T(counter1) = T_num/T_denum;
        counter1=counter1+1;

        if mod(jj-angle(z_low), 2*pi)==0
            index(counter2) = -((jj/pi-1)*((length(phi)-1)/((phi(1)-phi(end))/pi))-1);
            counter2 = counter2+1;
        end
    end

    func = c2*mean(T((index(1)-BW/2):(index(1)+BW/2))) + c2*mean(T((index(2)-BW/2):(index(2)+BW/2))) - c1*mean(T(1:(index(1)-BW/2)))- c1*mean(T((index(1)+BW/2):(index(2)-BW/2)))- c1*mean(T((index(2)+BW/2):end));
end