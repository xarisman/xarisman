clear all; %clc;
%d=16 tl=13 SLL = 20.52 th_st = 25
%d=14 tl=13 SLL = 18 th_st = 30-tol
%d=12 tl=12 SLL = 13.9
%d=10 tl = 15 SLL = 7.4
%d=8 tl = 11 SLL = 3.76
%d=6 tl = 16 SLL = 1.1

%tol = d-1
tic
% d = 12; N = 5; M=24;
% 
% th_1 = 30;
% th = [th_1, th_1+d, th_1+2*d, th_1+3*d, th_1+4*d, th_1+5*d];
% th_0 = th(2);
% th(2) = [];
% 
% th2 = [th_0, th]

%th3 = zeros(:, N+1);
p=0;
%theta2 = zeros(426, M);
for iiii=1:6
    d = 6+(iiii-1)*2;
    N = 5; M=24;
    th_1 = 30;
    counter = 0;
    for jj=1:1000
        for ii=1:N+1
            th = [th_1, th_1+d, th_1+2*d, th_1+3*d, th_1+4*d, th_1+5*d];
            th_0 = th(ii);
            th(ii) = [];
            th2 = [th_0, th];
            th3((N+1)*(jj-1)+ii, :) = th2;
        end
        
        th_1 = th_1+1;
        if max(th2)>=149
            break;
        end
        counter = counter+1;
    end
     
    A = zeros(M, M);
    Ai = zeros(M, M-1);
    ad = zeros(M, 1);
    e1 = zeros(1, M);
    e1(1) = 1;
    e1 = transpose(e1);
    a_th = zeros(M, 1);
    counter2 = 0;
    counter4=0;
    counter3 = 0;
    Rgg = eye(M);
    Rgigi = eye(M-1);
    s_2 = 0.01;
    Rnn = s_2*eye(M);
    
    %th_nsb = theta2(105, :);
    
    AF = zeros(1, 1800);
    th = 1:1800;
    
    I = eye(M);
    SINR_db = zeros(1, length(th3(:, :)));
    SLL = zeros(1, length(th3(:, :)));
    dth = zeros(length(th3(:, :)), N+0);
    %disp(length(th3(:)))
    theta2 = zeros(length(th3(:, :)), M);
    for n=1:length(th3(:, :))
        
        p=0;
        theta = th3(n, :);
        theta_0 = theta(1);
        %d6d8d10d12d14-t11-1.3,4.32,8.8,13.7,17.46 d16-t13-20.78
        if d>13
            tol = 13;
        elseif d==6
            tol =21;
        else
            tol =15;
        end
        tol2 = d/100;
        if d==6
            tol2 = d/20;
        end

        
        
        if theta_0<=75
            th_st = 37-tol+0.5; 
        elseif theta_0>125
            th_st = 33-tol+0.5;
        else
            th_st = 35-tol+0.5;
        end
        total_range = 180-2.05*tol-2*th_st-1.4*tol2;
        
        if d~=6
            num_parts = 18;
            num_parts_1 = 9;
            num_parts_2 = num_parts-num_parts_1;
            part_length = total_range/num_parts;
            angles = zeros(1, num_parts);
        elseif d==8
            num_parts = 14;
            num_parts_1 = 7;
            num_parts_2 = num_parts-num_parts_1;
            part_length = total_range/num_parts;
            angles = zeros(1, num_parts);
        elseif d==12
            num_parts = 15;
            num_parts_1 = 9;
            num_parts_2 = num_parts-num_parts_1;
            part_length = total_range/num_parts;
            angles = zeros(1, num_parts);
        else
            num_parts = 15;
            num_parts_1 = 8;
            num_parts_2 = num_parts-num_parts_1;
            part_length = total_range/num_parts;
            angles = zeros(1, num_parts);
        end

        if theta_0<=60
            total_range1 = total_range/2+5;
            total_range2 = total_range-total_range1;
        elseif theta_0>130
            total_range1 = total_range/2-3;
            total_range2 = total_range-total_range1;
        else
            total_range1 = total_range/2;
            total_range2 = total_range-total_range1;
        end

        part_length1 = total_range1/num_parts_1;
        part_length2 = total_range2/num_parts_2;

        tr=0;
        cc=1;
        angles(1) = 14;
        angles(18) = 156;
        if theta_0<42
            angles(1) = 1;
            angles(18) = 156;
        end
        if theta_0>138
            angles(1) = 28;
            angles(18) = 170;
        end

        for i = 2:num_parts-1
            if i <total_range1
                part_length = part_length1;
            else
                part_length = part_length2;
            end
            angles(i) = th_st + (-1+i)*part_length+2*tol*tr;
            %angles(i) = (i+1)*part_length+2*tol*tr;
            if abs(angles(i)-theta_0)<1*tol
                angles(i) = angles(i)+2.05*tol;
                tr=1;
            elseif abs(angles(i)-theta_0)<40
                %part_length = part_length1-2;
            end
        
            if abs(angles(i)-th3(n, cc+1))<tol2
                angles(i) = angles(i)+1*tol2;
                cc = cc+1;
            end
            
        end
        theta_0 = theta(1);
        %disp(theta)
        if d==6
            for ii=1:N
                if theta(ii)<theta_0
                    theta(ii)=theta(ii)-tol;
                elseif theta(ii)>theta_0
                    theta(ii)= theta(ii)+tol;
                else
                    theta(ii)=theta(ii);
                end
            end
            theta = [theta, angles];
        elseif d==8
            for ii=1:N
                if theta(ii)<theta_0
                    theta(ii)=theta(ii)-1;
                elseif theta(ii)>theta_0
                    theta(ii)= theta(ii)+1;
                else
                    theta(ii)=theta(ii);
                end
            end
            theta = [theta, angles];
        elseif d==10
            for ii=1:N
                if abs(theta(ii)-theta_0)<10 && theta(ii)<theta_0
                    theta(ii)=theta(ii)-1;
                elseif abs(theta(ii)-theta_0)<10 && theta(ii)>theta_0
                    theta(ii)= theta(ii)+1;
                else
                    theta(ii)=theta(ii);
                end
            end
            theta = [theta, angles];
        else
            theta = [theta, angles];
        end
        %theta_0 = theta(1);
        theta = sort(theta(2:end));
        %theta = theta(2:end);
        theta = [theta_0, theta];
        %disp(theta)
        theta2(n, :) = theta;
    
        counter3 = counter3+1;
        th_nsb = theta2(n, :);
        %disp(th_nsb(1))
        for ii=1:M
            for jj=1:M
                A(ii, jj) = exp((1j)*(ii-1)*pi*cosd(th_nsb(jj)));
                if jj>1
                    Ai(ii, jj-1) = exp((1j)*(ii-1)*pi*cosd(th_nsb(jj)));
                end
            end
        end
        for ii=1:M
            ad(ii) = exp((1i)*pi*ii*cosd(th_nsb(1)));
        end
        A_H = transpose(A);
        A_Hi = transpose(Ai);
        w_nsb = (A_H+s_2*I)\e1;
        
        Rxx = A*Rgg*A_H + Rnn;
        Ruu = Ai*Rgigi*A_Hi + Rnn;
        
        l_max = abs(eigs(Ruu\Rxx, 1));
%         SINR = l_max+1;
%         SINR_db(n) = 10*log10(SINR);
        SINR = 1*transpose(w_nsb)*ad*transpose(ad)*w_nsb/(transpose(w_nsb)*Ai*Rgigi*A_Hi*w_nsb+transpose(w_nsb)*Rnn*w_nsb);
        SINR_db(n) = 10*log10(abs(SINR));
        
        for th=1:1800
            for ii=1:M
                a_th(ii) = exp((1j)*pi*(ii-1)*cosd(th/10));
            end
            AF(th) = transpose(w_nsb)*a_th;
        end
        x = 1:length(AF);
        AF = AF/max(abs(AF));
    
        if(n<180&& n>150 && d==6 && 1==1)
            figure
            plot(x/10, abs(AF))
            hold on
            xline(th_nsb(1), 'r')
            hold off
            
        end
    
        %nullIndices = find(AF == 0);
        
        [null, delta_th_0] = max(AF);
        if abs(delta_th_0/10-th_nsb(1))>tol
            p=1;

            counter2 = counter2+1;
        end
        counter4=counter4+1;
        
        %SLL
        sll_array = [AF(1:(th_nsb(1)*10-10*d)) AF((th_nsb(1)*10+10*d):end)];
        [maxValue, maxIndex] = max(AF);
        second_max_value = abs(max(sll_array));
    
        SLL(n) = 20*log10(abs(maxValue)/abs(second_max_value));
        
        
        diff_A = diff(AF);
        local_min_indices = find(diff_A(1:end-1) > 0 & diff_A(2:end) < 0) + 1;
        
        
        for x=1:N
            [~, dth(n, x)] = min(abs(local_min_indices-10*th3(n, x+1)));
            dth(n, x) = dth(n, x)/10;
        end
        
        [~, index_max] = max(AF(n));
        dth(6,index_max) = abs(index_max-10*th3(1))/10;
    end
    %disp(SLL(1))
    disp("d: " + d)
    disp("mean(SLL): " + mean(SLL))
    disp("SINR_db_mean: " + mean(SINR_db))
    disp("counter2: " + counter2 + "/" + counter4 +", "+ counter2/counter4*100 +"%")
    disp("th0: " + mean(dth(6,:)))
    disp("dth1: " + mean(dth(1,:)))
    disp("dth2: " + mean(dth(2,:)))
    disp("dth3: " + mean(dth(3,:)))
    disp("dth4: " + mean(dth(4,:)))
    disp("dth5: " + mean(dth(5,:)))
    
    if iiii<6
        clear all;
    end
end
toc