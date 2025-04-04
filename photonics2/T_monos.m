function [Tm] = T_monos(twi, ai, start_phi, end_phi, num_points)
    tw = twi;
    a = ai;
    phi = linspace(start_phi, end_phi, num_points);
    counter1=1;
    T = zeros(1, length(phi));
    
    
    for i=phi
        T(1, counter1) = (tw^2+a^2-2*a*tw*cos(i))/(1+(a*tw)^2-2*a*tw*cos(i));
        counter1 = counter1+1;
    end
    
    Tm = T(:);
end