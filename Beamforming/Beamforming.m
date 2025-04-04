
delta=6;
SNR=20;

initial_angle = 30;
end_angle= 150;
max1 = end_angle-5*delta;

Data = zeros(6*(max1-initial_angle+1),14);
j=1;

while(initial_angle<=max1)
    
    
    theta_array=six_sets(initial_angle,delta);
    
    for i=1:6
        
        theta = theta_array(i,:);
        
        Data(j,1)=theta(1);
        Data(j,2)=theta(2);
        Data(j,3)=theta(3);
        Data(j,4)=theta(4);
        Data(j,5)=theta(5);
        Data(j,6)=theta(6);
        
        result = NSB_with_angles(SNR,theta);
        
        Data(j,7)=result(1);
        Data(j,8)=result(2);
        Data(j,9)=result(3);
        Data(j,10)=result(4);
        Data(j,11)=result(5);
        Data(j,12)=result(6);
        
        Data(j,13)=result(7);
        Data(j,14)=result(8);
        
        j=j+1;
    end
    
    initial_angle=initial_angle+1;
end


fileID=fopen('AoAdev_SINR_SLL.txt','w');
formatSpec=('%.3f ,%.3f ,%.3f ,%.3f ,%.3f, %.3f , %.3f , %.3f , %.3f , %.3f , %.3f , %.3f , %.3f , %.3f \n');
fprintf(fileID,formatSpec,Data');
fclose(fileID);


%main lobe divergence:
min_MLB=min(Data(:,7))
max_MLB=max(Data(:,7))
mean_MLB=mean(Data(:,7))
std_MLB=std(Data(:,7))

%null divergence:
min_ND=min([Data(:,8),Data(:,9),Data(:,10),Data(:,11),Data(:,12)])
max_ND=max([Data(:,8),Data(:,9),Data(:,10),Data(:,11),Data(:,12)])
mean_ND=mean([Data(:,8),Data(:,9),Data(:,10),Data(:,11),Data(:,12)])
std_ND=std([Data(:,8),Data(:,9),Data(:,10),Data(:,11),Data(:,12)])

%SINR:
min_SINR=min(Data(:,13))
max_SINR=max(Data(:,13))
mean_SINR=mean(Data(:,13))
std_SINR=std(Data(:,13))

%SLL:
min_SLL=min(Data(:,14))
max_SLL=max(Data(:,14))
mean_SLL=mean(Data(:,14))
std_SLL=std(Data(:,14))



    
