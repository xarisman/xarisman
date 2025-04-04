function [final_array] = adding_angles(original_array)
    %this function is used to add the fictitious angles

    mainlobe = original_array(1);  
    hpbw = 4; 

    num_fake_angles = 18; 
    num_total = numel(original_array) + num_fake_angles; 

    lower_bound = 0; 
    upper_bound = 180;

    spacing = (upper_bound - lower_bound) / (num_total - 1);  

    range1 = linspace(lower_bound, mainlobe - hpbw, ceil(num_fake_angles/2));
    range2 = linspace(mainlobe + hpbw, upper_bound, floor(num_fake_angles/2));
    fictitious_angles = [range1, range2];

    final_array = [original_array, fictitious_angles];  

end