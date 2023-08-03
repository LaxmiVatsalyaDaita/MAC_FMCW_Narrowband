%AFter back-off period with prob p the radar enters into sensing state 
%It senses the channel for f_min, if it finds any it will enter into
%back-off in the next slot or start transmittiong chirp.
%the states are 0(wait), -1(CCA1),-2(CCA2), 1,...K, transmit
%fc1 = f_max - \Deltaf, fc2 = f_min

function theta = CSMA_noncollocated_twocca_real(N,K,L,p,T)
CCA1_fail = zeros(1,N);
CCA2_fail = zeros(1,N);
coll_ind = zeros(N,N);
coll =  zeros(1,N);
packet_cycle = zeros(1,N);
sense1_flag =  zeros(1,N);
sense2_flag = zeros(1,N);
sensed1 = zeros(N);
sensed2 = zeros(N);
D = 0;
X = zeros(D+2,N);
X(1,:) = randi([1 K*L],1,N);
del = rand(1,N); %clock offsets with respect to an universal clock
[~, sorted_id] = sort(del);
delta = repmat(del',1,N) - repmat(del,N,1); %relative clock offset matrix
u = D*triu(rand(N));
d = u + u'- 2*diag(u(1:N+1:N^2)); 
tau = floor(d + delta);
%Random states correspond to initial dummy packet transmissions. 
%Actual activities start at random times, at the end of dummy packets.%del = [0.9842    0.6284    0.5694    0.8980    0.2170    0.4578    0.6538    0.1917    0.7327    0.9100]

coll_time = K*L-X(1,:); %This is to avoid overcounting of collisions
for t = 1:T
    X(2:D+2,:) = X(1:D+1,:);
    for i = 1:N
        if (X(2,i) == K*L)||(X(2,i) == 0)
            X(1,i) = -1*binornd(1,p);
        end
        if (X(2,i) > 0)&&(X(2,i) < K*L)
            X(1,i) = X(2,i)+1;
        end
        if (X(2,i) == -1) && (sense1_flag(i) == 0)
            X(1,i) = -2;
        end
        if (X(2,i) == -1) && (sense1_flag(i) == 1)
            X(1,i) = -1*binornd(1,p);                        
            CCA1_fail(i) = CCA1_fail(i)+1;
        end        
        if (X(2,i) == -2) && (sense2_flag(i) == 0)
            X(1,i) = 1;
            packet_cycle(i) = packet_cycle(i) + 1;    
        end
        if (X(2,i) == -2) && (sense2_flag(i) == 1)
            X(1,i) = -1*binornd(1,p);                        
            CCA2_fail(i) = CCA2_fail(i)+1;
        end        
    end
    % X
    % packet_cycle
    
    if t < K*L
        for n = 1:N
            i = sorted_id(n);
            for j = [1:i-1 i+1:N]
                if tau(j,i) >= 0
                    sensed1(i,j) = (X(tau(j,i)+1,j) > 0)*(X(tau(j,i)+1,j) <= (K-1)*L)*((mod(X(tau(j,i)+1,j),K)==K-1)+(mod(X(tau(j,i)+1,j),K) == 0))*(packet_cycle(j) > 0);
                else 
                    sensed1(i,j) = (X(1,j) > 0)*(X(1,j) <= (K-1)*L)*((mod(X(1,j),K)==K-2)+(mod(X(1,j),K) == K-1))*(packet_cycle(j) > 0);
                end
                if tau(j,i) >= 0
                    sensed2(i,j) = (mod(X(tau(j,i)+1,j),K)==1)*(packet_cycle(j) > 0);
                else 
                    sensed2(i,j) = (X(1,j) == -2)*(sense2_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0)*(packet_cycle(j) > 0);
                end
            end
            sense1_flag(i) = sum(sensed1(i,1:N)) > 0;
            sense2_flag(i) = sum(sensed2(i,1:N)) > 0;
        end
    else
        for n = 1:N
            i = sorted_id(n);
            for j = [1:i-1 i+1:N]
                if tau(j,i) >= 0
                    sensed1(i,j) = (X(tau(j,i)+1,j) > 0)*(X(tau(j,i)+1,j) <= (K-1)*L)*((mod(X(tau(j,i)+1,j),K)==K-1)+(mod(X(tau(j,i)+1,j),K) == 0));
                else 
                    sensed1(i,j) = (X(1,j) > 0)*(X(1,j) <= (K-1)*L)*(mod(X(1,j),K)==K-2)+(mod(X(1,j),K) == K-1);
                end
                if tau(j,i) >= 0
                    sensed2(i,j) = mod(X(tau(j,i)+1,j),K)==1;
                else 
                    sensed2(i,j) = (X(1,j) == -2)*(sense2_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0);
                end
            end
            sense1_flag(i) = sum(sensed1(i,1:N)) > 0;
            sense2_flag(i) = sum(sensed2(i,1:N)) > 0;
        end
    end
%    sense1_flag
%    sense2_flag

    if t < K*L
        for n = 1:N
            i = sorted_id(n);
            if (mod(X(1,i),K) == 1)&&(t > coll_time(i)-1)  
                for j = [1:i-1 i+1:N]
                    if tau(j,i) >= 0
                        coll_ind(j,i) = ((mod(X(tau(j,i)+1,j),K) == 1) + (mod(X(tau(j,i)+2,j),K) == 1))*(packet_cycle(j) > 0);    
                    else
                        coll_ind(j,i) = (X(1,j) == -2)*(sense2_flag(j) == 0) ... 
                            + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0)*(packet_cycle(j) > 0) + (mod(X(1,j),K) == 1)*(packet_cycle(j) > 0);
                    end
                end
                if sum(coll_ind(:,i)) > 0
                    coll(i) = coll(i) + 1;
                    coll_time(i) = t+K*L-X(1,i);
                end
            end
        end 
    else
        for n = 1:N
            i = sorted_id(n);
            if (mod(X(1,i),K) == 1)&&(t > coll_time(i)-1)  
                for j = [1:i-1 i+1:N]
                    if tau(j,i) >= 0
                        coll_ind(j,i) = (mod(X(tau(j,i)+1,j),K) == 1) + (mod(X(tau(j,i)+2,j),K) == 1);    
                    else
                        coll_ind(j,i) = (X(1,j) == -2)*(sense2_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0) + (mod(X(1,j),K) == 1);
                    end
                end
                if sum(coll_ind(:,i)) > 0
                    coll(i) = coll(i) + 1;
                    coll_time(i) = t+K*L-X(1,i);
                end
            end
        end
    end
%  coll
end 
CCA1_fail;
CCA2_fail;
packet_cycle;
coll;
% mean_CCA_fail = mean(CCA_fail);
% mean_packet_cycles = mean(packet_cycle);
% mean_coll = mean(coll);
theta = (sum(packet_cycle)-sum(coll))*K*L/T
end