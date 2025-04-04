function [combos]= six_sets(initial_angle,delta)
    %creates all the six sets of angles starting from the initial angle and sets them apart by delta
    
    combos(1,:) = [initial_angle, initial_angle+delta,initial_angle+2*delta,initial_angle+3*delta,initial_angle+4*delta,initial_angle+5*delta];
    
    for i=1:5

        combos(i+1,:) = combos(i,:);
        [combos(i+1,1),combos(i+1,i+1)] = deal(combos(i+1,i+1),combos(i+1,1));

    end


end