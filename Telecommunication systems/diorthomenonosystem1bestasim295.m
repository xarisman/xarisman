N = 10^5; %number of symbols to process
M = 4;
N2 = 100;
SNR = -3:20;
ipHat = zeros(1,N); %creates array with N datatype
ipHat2 = zeros(1,N);
array_3 = zeros(1,100);
array_1 = zeros(1,length(SNR));

for ii = 1:length(SNR)
    for i = 1:100
        a = (i-1)./100;
        alpha4pam = [-3*sqrt(a./5) -sqrt(a./5) sqrt(a./5) 3*sqrt(a./5)]; %cordinates of 4-PAM
        alpha4pam2 = [-3*sqrt((1-a)./5) -sqrt((1-a)./5) sqrt((1-a)./5) 3*sqrt((1-a)./5)];
        
        %Simulation of U1
        ip = randsrc(1, N, alpha4pam);
        s = (1/sqrt(a)).*ip; %normaliazation of energy
        n = 1/sqrt(a.*4).*(randn(1,N) + (1j)*randn(1,N));
        
        y = s + 10^(-SNR(ii)./20).*n; % additive white gaussian noise
    
        %demodulation
        r = real(y);
        %areas of decisions
        ipHat(find(r< -2/sqrt(5))) = -3*sqrt(a./5);
        ipHat(find(r>= 2/sqrt(5))) = 3*sqrt(a./5);
        ipHat(find(r>=-2/sqrt(5) & r<0)) = -sqrt(a./5);
        ipHat(find(r>=0 & r<2/sqrt(5))) = sqrt(a/5);
        
        nErr = size(find([ip - ipHat]), 2); %counting errors
        
        %Simulation of U2
        ip2 = randsrc(1, N, alpha4pam2);
        s2 = (1/sqrt(1-a)).*ip2; %normaliazation of energy
        n2 = 1/sqrt((1-a).*8)*(randn(1,N) + (1j)*randn(1,N));
        
        y2 = s2 + 10^(-SNR(ii)./20).*n2; % additive white gaussian noise
    
        %demodulation
        r2 = real(y2);
        %areas of decisions
        ipHat2(find(r2< -2/sqrt(5))) = -3*sqrt((1-a)./5);
        ipHat2(find(r2>= 2/sqrt(5))) = 3*sqrt((1-a)./5);
        ipHat2(find(r2>=-2/sqrt(5) & r2<0)) = -sqrt((1-a)./5);
        ipHat2(find(r2>=0 & r2<2/sqrt(5))) = sqrt((1-a)./5);
        
        nErr2 = size(find([ip2 - ipHat2]), 2); %counting errors

        array_3(i) = max(nErr/N, nErr2/N); %Saving the Maximum between SEP1 and SEP2
       
    end
    [M, I] = min(array_3); %Saving the minimum value of error in every SNR and then saving the index
    array_1(ii) = I;
    
end

plot(array_1)