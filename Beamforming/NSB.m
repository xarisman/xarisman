function [result] = NSB(SNR, theta)
    %this is the null steering beamforming algorithm as given without
    %adding fictitious angles 
    
    theta_Rad = theta*pi/180;
    
    P=1;
    Pn=10^(-SNR/10);
    
    M=24;
    N=length(theta)-1;
    
    Rnn = Pn*eye(M);
    Rgg=P*eye(N+1);
    
    e1=zeros(N+1,1);
    e1(1,1)=1;
    
    A=zeros(M,N+1);
    
    for i=1:M
        for j=1:N+1
            A(i,j)=exp(1i*pi*(i-1)*cos(theta_Rad(j)));
        end 
    end
    
    a_n = zeros(M,1);
    for i=1:M
        a_n(i,1)=exp(1i*pi*(i-1)*cos(theta_Rad(1)));
    end
    
    
    Rxx=A*Rgg*A'+Rnn;
    
    
    
    Pn=0.001;
    
    determinant=abs(det(A'*A));
    
    if abs(det(A'*A))<=10^(-6)
        W_NSB=A*inv(A'*A+Pn*eye(N+1))*e1;
    else
        W_NSB=A*inv(A'*A)*e1;
    end
    
    
    
    theta180=0:0.1:180;
    theta180_Rad=theta180*pi/180;
    
    
    a_theta=zeros(M,length(theta180_Rad));
    for i=1:M
        a_theta(i,:)=exp(1i*pi*(i-1)*cos(theta180_Rad));
    end
    
    AF=W_NSB'*a_theta;
    
    figure(1)
    plot(theta180,abs(AF))
    xlabel('Elevation Angle')
    ylabel('Magnitute of AF')
    findpeaks(abs(AF))
    
    [peaks, locations] = findpeaks(abs(AF));
    
    locations2 = locations*0.1;
    
    deviations = zeros(N+1,1);
    
    for i=1:N+1
        
        deviations(i)=min(abs(locations2-theta(i)));
        
    end
    
    %SINR=real((W_NSB'*a_n*a_n'*W_NSB)/(W_NSB'*A*Rgg*A'*W_NSB));
    
    Rgigi=eye(5,5);
    
    Ai=A(:,2:N+1);
    Ruu=Ai*Rgigi*Ai'+Rnn;
    %SINR = real((W_NSB'*Rxx*W_NSB)/(W_NSB'*Ruu*W_NSB)-1);
    
    SINR=real((W_NSB'*a_n*a_n'*W_NSB)/(W_NSB'*Ai*Rgigi*Ai'*W_NSB+W_NSB'*Rnn*W_NSB));
    SINR_dB= 10*log10(SINR);
    
    biggestpeaks = maxk(peaks,2);
    SLL = min(biggestpeaks);
    
    SLL_dB = 10*log10(SLL);
    
    result=zeros(N+3,1);
    result(1)=deviations(1);
    result(2)=deviations(2);
    result(3)=deviations(3);
    result(4)=deviations(4);
    result(5)=deviations(5);
    result(6)=deviations(6);
    result(7)=SINR_dB;
    result(8)=SLL_dB;
    
end