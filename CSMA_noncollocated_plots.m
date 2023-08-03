%AFter back-off period with prob p the radar enters into sensing state 
%It senses the channel for f_min, if it finds any it will enter into
%back-off in the next slot or start transmittiong chirp.
%the states are 0(wait), -1(sense), 1,...K, transmit

%assumption that the radars are not synchronized with each other but are
%are all in sync with a universal clock. All the ith slots fall within
%delta.

%function theta = CSMA_noncollocated_plots(N,K,L,p,T,D) %N is no. of nodes, K is no. of slots per chirp, L is no. of chirps per packet, p is probability with which CCA is done in the (K*L+1)th slot 
the=[];
N=50;
K=40;
T=10000;
D=0;
L=32;
p_values = [0:0.1:1];

X_(1,:) = randi([1 K*L],1,N);
del_ = rand(1,N); %clock offsets with respect to an universal clock
%u = D*triu(rand(N));
tri = triu(rand(N));



for xyz=1:length(p_values)
    
    CCA_fail = zeros(1,N);
    coll_ind = zeros(N,N);
    coll =  zeros(1,N);
    packet_cycle = zeros(1,N); % for counting number of packet transmissions
    sense_flag =  zeros(1,N);
    sensed = zeros(N);
    %D = 100; % distance between nodes.
    X = zeros(D+2,N);
    %X(1,:) = randi([1 K*L],1,N);
    X(1,:) = X_(1,:);
    del = del_; %clock offsets with respect to an universal clock
    [~, sorted_id] = sort(del); % sorting the clock offsets in increasing order
    %del = sort(del);
    delta = repmat(del',1,N) - repmat(del,N,1); %relative clock offset matrix
    %u = D*triu(rand(N));
    u = D*tri;
    d = u + u'- 2*diag(u(1:N+1:N^2)); % symmetric matrix with diagonal elements = 0
    tau = floor(d + delta);
    %Random states correspond to initial dummy packet transmissions.
    %Actual activities start at random times, at the end of dummy packets.
    coll_time = K*L-X(1,:); %This is avoid overcounting of collisions
    for t = 1:T
        X(2:D+2,:) = X(1:D+1,:);
        for i = 1:N
            if (X(2,i) == K*L)||(X(2,i) == 0) % last slot or backoff state
                X(1,i) = -1*binornd(1,p_values(xyz)); % CCA with probability p
            end
            if (X(2,i) > 0)&&(X(2,i) < K*L) % packet transmission states i.e 1 to K*L
                X(1,i) = X(2,i)+1; % incerement to next state
            end
            if (X(2,i) == -1) && (sense_flag(i) == 0) % CCA success
                X(1,i) = 1; % start transmission
                packet_cycle(i) = packet_cycle(i) + 1; % increment packets transmitted
            end
            if (X(2,i) == -1) && (sense_flag(i) == 1) % CCA fail
                X(1,i) = -1*binornd(1,p_values(xyz)); % do CCA again with probability p
                CCA_fail(i) = CCA_fail(i)+1; % increment the count of CCA fails
            end
        end
        % X
        % packet_cycle

        if t < K*L % the packet being transmitted may be actual or virtual packet; have to check the packet cycle count
            for n = 1:N % for all nodes
                i = sorted_id(n);
                % i and j are nodes
                for j = [1:i-1 i+1:N] % prev node and next node
                    if tau(j,i) >= 0 %
                        sensed(i,j) = (mod(X(tau(j,i)+1,j),K)==1)*(packet_cycle(j) > 0);
                    else
                        sensed(i,j) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0)*(packet_cycle(j) > 0);
                    end
                end
                sense_flag(i) = sum(sensed(i,1:N)) > 0;
            end
            %        sensed
        else
            for n = 1:N
                i = sorted_id(n);
                for j = [1:i-1 i+1:N]
                    if tau(j,i) >= 0
                        sensed(i,j) = mod(X(tau(j,i)+1,j),K)==1;
                    else
                        sensed(i,j) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0);
                    end
                end
                sense_flag(i) = sum(sensed(i,1:N)) > 0;
            end
        end
        %    sense_flag;

        if t < K*L
            for n = 1:N
                i = sorted_id(n);
                if (mod(X(1,i),K) == 1)&&(t > coll_time(i)-1)
                    for j = [1:i-1 i+1:N]
                        if tau(j,i) >= 0
                            coll_ind(j,i) = ((mod(X(tau(j,i)+1,j),K) == 1) + (mod(X(tau(j,i)+2,j),K) == 1))*(packet_cycle(j) > 0);
                        else
                            coll_ind(j,i) = (X(1,j) == -1)*(sense_flag(j) == 0) ...
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
                            coll_ind(j,i) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0) + (mod(X(1,j),K) == 1);
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
    CCA_fail;
    packet_cycle;
    coll;
    % mean_CCA_fail = mean(CCA_fail);
    % mean_packet_cycles = mean(packet_cycle);
    % mean_coll = mean(coll);

    theta = (sum(packet_cycle)-sum(coll))*K*L/T
    the(end+1)=theta;
    %print(p_values(j));

    %plot(p_values,theta)
    %end
end
%hold off
the;
figure;
%subplot(2,3,1);
plot(p_values,the,'r')
%title("L=" + L + " , D=0")

hold on

the=[];
D=1;

for xyz=1:length(p_values)
    
    CCA_fail = zeros(1,N);
    coll_ind = zeros(N,N);
    coll =  zeros(1,N);
    packet_cycle = zeros(1,N); % for counting number of packet transmissions
    sense_flag =  zeros(1,N);
    sensed = zeros(N);
    %D = 100; % distance between nodes.
    X = zeros(D+2,N);
    %X(1,:) = randi([1 K*L],1,N);
    X(1,:) = X_(1,:);
    del = del_; %clock offsets with respect to an universal clock
    [~, sorted_id] = sort(del); % sorting the clock offsets in increasing order
    %del = sort(del);
    delta = repmat(del',1,N) - repmat(del,N,1); %relative clock offset matrix
    %u = D*triu(rand(N));
    u = D*tri;
    d = u + u'- 2*diag(u(1:N+1:N^2)); % symmetric matrix with diagonal elements = 0
    tau = floor(d + delta);
    %Random states correspond to initial dummy packet transmissions.
    %Actual activities start at random times, at the end of dummy packets.
    coll_time = K*L-X(1,:); %This is avoid overcounting of collisions
    for t = 1:T
        X(2:D+2,:) = X(1:D+1,:);
        for i = 1:N
            if (X(2,i) == K*L)||(X(2,i) == 0) % last slot or backoff state
                X(1,i) = -1*binornd(1,p_values(xyz)); % CCA with probability p
            end
            if (X(2,i) > 0)&&(X(2,i) < K*L) % packet transmission states i.e 1 to K*L
                X(1,i) = X(2,i)+1; % incerement to next state
            end
            if (X(2,i) == -1) && (sense_flag(i) == 0) % CCA success
                X(1,i) = 1; % start transmission
                packet_cycle(i) = packet_cycle(i) + 1; % increment packets transmitted
            end
            if (X(2,i) == -1) && (sense_flag(i) == 1) % CCA fail
                X(1,i) = -1*binornd(1,p_values(xyz)); % do CCA again with probability p
                CCA_fail(i) = CCA_fail(i)+1; % increment the count of CCA fails
            end
        end
        % X
        % packet_cycle

        if t < K*L % the packet being transmitted may be actual or virtual packet; have to check the packet cycle count
            for n = 1:N % for all nodes
                i = sorted_id(n);
                % i and j are nodes
                for j = [1:i-1 i+1:N] % prev node and next node
                    if tau(j,i) >= 0 %
                        sensed(i,j) = (mod(X(tau(j,i)+1,j),K)==1)*(packet_cycle(j) > 0);
                    else
                        sensed(i,j) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0)*(packet_cycle(j) > 0);
                    end
                end
                sense_flag(i) = sum(sensed(i,1:N)) > 0;
            end
            %        sensed
        else
            for n = 1:N
                i = sorted_id(n);
                for j = [1:i-1 i+1:N]
                    if tau(j,i) >= 0
                        sensed(i,j) = mod(X(tau(j,i)+1,j),K)==1;
                    else
                        sensed(i,j) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0);
                    end
                end
                sense_flag(i) = sum(sensed(i,1:N)) > 0;
            end
        end
        %    sense_flag;

        if t < K*L
            for n = 1:N
                i = sorted_id(n);
                if (mod(X(1,i),K) == 1)&&(t > coll_time(i)-1)
                    for j = [1:i-1 i+1:N]
                        if tau(j,i) >= 0
                            coll_ind(j,i) = ((mod(X(tau(j,i)+1,j),K) == 1) + (mod(X(tau(j,i)+2,j),K) == 1))*(packet_cycle(j) > 0);
                        else
                            coll_ind(j,i) = (X(1,j) == -1)*(sense_flag(j) == 0) ...
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
                            coll_ind(j,i) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0) + (mod(X(1,j),K) == 1);
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
    CCA_fail;
    packet_cycle;
    coll;
    % mean_CCA_fail = mean(CCA_fail);
    % mean_packet_cycles = mean(packet_cycle);
    % mean_coll = mean(coll);

    theta = (sum(packet_cycle)-sum(coll))*K*L/T;
    the(end+1)=theta;
    %print(p_values(j));

    %plot(p_values,theta)
    %end
end
%hold off
the;
%subplot(2,3,2);
plot(p_values,the,'c')
%title("L=" + L + " , D=1")




the=[];
D=2;

for xyz=1:length(p_values)
    
    CCA_fail = zeros(1,N);
    coll_ind = zeros(N,N);
    coll =  zeros(1,N);
    packet_cycle = zeros(1,N); % for counting number of packet transmissions
    sense_flag =  zeros(1,N);
    sensed = zeros(N);
    %D = 100; % distance between nodes.
    X = zeros(D+2,N);
    %X(1,:) = randi([1 K*L],1,N);
    X(1,:) = X_(1,:);
    del = del_; %clock offsets with respect to an universal clock
    [~, sorted_id] = sort(del); % sorting the clock offsets in increasing order
    %del = sort(del);
    delta = repmat(del',1,N) - repmat(del,N,1); %relative clock offset matrix
    %u = D*triu(rand(N));
    u = D*tri;
    d = u + u'- 2*diag(u(1:N+1:N^2)); % symmetric matrix with diagonal elements = 0
    tau = floor(d + delta);
    %Random states correspond to initial dummy packet transmissions.
    %Actual activities start at random times, at the end of dummy packets.
    coll_time = K*L-X(1,:); %This is avoid overcounting of collisions
    for t = 1:T
        X(2:D+2,:) = X(1:D+1,:);
        for i = 1:N
            if (X(2,i) == K*L)||(X(2,i) == 0) % last slot or backoff state
                X(1,i) = -1*binornd(1,p_values(xyz)); % CCA with probability p
            end
            if (X(2,i) > 0)&&(X(2,i) < K*L) % packet transmission states i.e 1 to K*L
                X(1,i) = X(2,i)+1; % incerement to next state
            end
            if (X(2,i) == -1) && (sense_flag(i) == 0) % CCA success
                X(1,i) = 1; % start transmission
                packet_cycle(i) = packet_cycle(i) + 1; % increment packets transmitted
            end
            if (X(2,i) == -1) && (sense_flag(i) == 1) % CCA fail
                X(1,i) = -1*binornd(1,p_values(xyz)); % do CCA again with probability p
                CCA_fail(i) = CCA_fail(i)+1; % increment the count of CCA fails
            end
        end
        % X
        % packet_cycle

        if t < K*L % the packet being transmitted may be actual or virtual packet; have to check the packet cycle count
            for n = 1:N % for all nodes
                i = sorted_id(n);
                % i and j are nodes
                for j = [1:i-1 i+1:N] % prev node and next node
                    if tau(j,i) >= 0 %
                        sensed(i,j) = (mod(X(tau(j,i)+1,j),K)==1)*(packet_cycle(j) > 0);
                    else
                        sensed(i,j) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0)*(packet_cycle(j) > 0);
                    end
                end
                sense_flag(i) = sum(sensed(i,1:N)) > 0;
            end
            %        sensed
        else
            for n = 1:N
                i = sorted_id(n);
                for j = [1:i-1 i+1:N]
                    if tau(j,i) >= 0
                        sensed(i,j) = mod(X(tau(j,i)+1,j),K)==1;
                    else
                        sensed(i,j) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0);
                    end
                end
                sense_flag(i) = sum(sensed(i,1:N)) > 0;
            end
        end
        %    sense_flag;

        if t < K*L
            for n = 1:N
                i = sorted_id(n);
                if (mod(X(1,i),K) == 1)&&(t > coll_time(i)-1)
                    for j = [1:i-1 i+1:N]
                        if tau(j,i) >= 0
                            coll_ind(j,i) = ((mod(X(tau(j,i)+1,j),K) == 1) + (mod(X(tau(j,i)+2,j),K) == 1))*(packet_cycle(j) > 0);
                        else
                            coll_ind(j,i) = (X(1,j) == -1)*(sense_flag(j) == 0) ...
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
                            coll_ind(j,i) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0) + (mod(X(1,j),K) == 1);
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
    CCA_fail;
    packet_cycle;
    coll;
    % mean_CCA_fail = mean(CCA_fail);
    % mean_packet_cycles = mean(packet_cycle);
    % mean_coll = mean(coll);

    theta = (sum(packet_cycle)-sum(coll))*K*L/T;
    the(end+1)=theta;
    %print(p_values(j));

    %plot(p_values,theta)
    %end
end
%hold off
the;
%subplot(2,3,3);
plot(p_values,the,'g')
%title("L=" + L + " , D=2")



the=[];
D=3;

for xyz=1:length(p_values)
    
    CCA_fail = zeros(1,N);
    coll_ind = zeros(N,N);
    coll =  zeros(1,N);
    packet_cycle = zeros(1,N); % for counting number of packet transmissions
    sense_flag =  zeros(1,N);
    sensed = zeros(N);
    %D = 100; % distance between nodes.
    X = zeros(D+2,N);
    %X(1,:) = randi([1 K*L],1,N);
    X(1,:) = X_(1,:);
    del = del_; %clock offsets with respect to an universal clock
    [~, sorted_id] = sort(del); % sorting the clock offsets in increasing order
    %del = sort(del);
    delta = repmat(del',1,N) - repmat(del,N,1); %relative clock offset matrix
    %u = D*triu(rand(N));
    u = D*tri;
    d = u + u'- 2*diag(u(1:N+1:N^2)); % symmetric matrix with diagonal elements = 0
    tau = floor(d + delta);
    %Random states correspond to initial dummy packet transmissions.
    %Actual activities start at random times, at the end of dummy packets.
    coll_time = K*L-X(1,:); %This is avoid overcounting of collisions
    for t = 1:T
        X(2:D+2,:) = X(1:D+1,:);
        for i = 1:N
            if (X(2,i) == K*L)||(X(2,i) == 0) % last slot or backoff state
                X(1,i) = -1*binornd(1,p_values(xyz)); % CCA with probability p
            end
            if (X(2,i) > 0)&&(X(2,i) < K*L) % packet transmission states i.e 1 to K*L
                X(1,i) = X(2,i)+1; % incerement to next state
            end
            if (X(2,i) == -1) && (sense_flag(i) == 0) % CCA success
                X(1,i) = 1; % start transmission
                packet_cycle(i) = packet_cycle(i) + 1; % increment packets transmitted
            end
            if (X(2,i) == -1) && (sense_flag(i) == 1) % CCA fail
                X(1,i) = -1*binornd(1,p_values(xyz)); % do CCA again with probability p
                CCA_fail(i) = CCA_fail(i)+1; % increment the count of CCA fails
            end
        end
        % X
        % packet_cycle

        if t < K*L % the packet being transmitted may be actual or virtual packet; have to check the packet cycle count
            for n = 1:N % for all nodes
                i = sorted_id(n);
                % i and j are nodes
                for j = [1:i-1 i+1:N] % prev node and next node
                    if tau(j,i) >= 0 %
                        sensed(i,j) = (mod(X(tau(j,i)+1,j),K)==1)*(packet_cycle(j) > 0);
                    else
                        sensed(i,j) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0)*(packet_cycle(j) > 0);
                    end
                end
                sense_flag(i) = sum(sensed(i,1:N)) > 0;
            end
            %        sensed
        else
            for n = 1:N
                i = sorted_id(n);
                for j = [1:i-1 i+1:N]
                    if tau(j,i) >= 0
                        sensed(i,j) = mod(X(tau(j,i)+1,j),K)==1;
                    else
                        sensed(i,j) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0);
                    end
                end
                sense_flag(i) = sum(sensed(i,1:N)) > 0;
            end
        end
        %    sense_flag;

        if t < K*L
            for n = 1:N
                i = sorted_id(n);
                if (mod(X(1,i),K) == 1)&&(t > coll_time(i)-1)
                    for j = [1:i-1 i+1:N]
                        if tau(j,i) >= 0
                            coll_ind(j,i) = ((mod(X(tau(j,i)+1,j),K) == 1) + (mod(X(tau(j,i)+2,j),K) == 1))*(packet_cycle(j) > 0);
                        else
                            coll_ind(j,i) = (X(1,j) == -1)*(sense_flag(j) == 0) ...
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
                            coll_ind(j,i) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0) + (mod(X(1,j),K) == 1);
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
    CCA_fail;
    packet_cycle;
    coll;
    % mean_CCA_fail = mean(CCA_fail);
    % mean_packet_cycles = mean(packet_cycle);
    % mean_coll = mean(coll);

    theta = (sum(packet_cycle)-sum(coll))*K*L/T;
    the(end+1)=theta;
    %print(p_values(j));

    %plot(p_values,theta)
    %end
end
%hold off
the;
%subplot(2,3,4);
plot(p_values,the,'y')
%title("L=" + L + " , D=3")



the=[];
D=4;

for xyz=1:length(p_values)
    
    CCA_fail = zeros(1,N);
    coll_ind = zeros(N,N);
    coll =  zeros(1,N);
    packet_cycle = zeros(1,N); % for counting number of packet transmissions
    sense_flag =  zeros(1,N);
    sensed = zeros(N);
    %D = 100; % distance between nodes.
    X = zeros(D+2,N);
    %X(1,:) = randi([1 K*L],1,N);
    X(1,:) = X_(1,:);
    del = del_; %clock offsets with respect to an universal clock
    [~, sorted_id] = sort(del); % sorting the clock offsets in increasing order
    %del = sort(del);
    delta = repmat(del',1,N) - repmat(del,N,1); %relative clock offset matrix
    %u = D*triu(rand(N));
    u = D*tri;
    d = u + u'- 2*diag(u(1:N+1:N^2)); % symmetric matrix with diagonal elements = 0
    tau = floor(d + delta);
    %Random states correspond to initial dummy packet transmissions.
    %Actual activities start at random times, at the end of dummy packets.
    coll_time = K*L-X(1,:); %This is avoid overcounting of collisions
    for t = 1:T
        X(2:D+2,:) = X(1:D+1,:);
        for i = 1:N
            if (X(2,i) == K*L)||(X(2,i) == 0) % last slot or backoff state
                X(1,i) = -1*binornd(1,p_values(xyz)); % CCA with probability p
            end
            if (X(2,i) > 0)&&(X(2,i) < K*L) % packet transmission states i.e 1 to K*L
                X(1,i) = X(2,i)+1; % incerement to next state
            end
            if (X(2,i) == -1) && (sense_flag(i) == 0) % CCA success
                X(1,i) = 1; % start transmission
                packet_cycle(i) = packet_cycle(i) + 1; % increment packets transmitted
            end
            if (X(2,i) == -1) && (sense_flag(i) == 1) % CCA fail
                X(1,i) = -1*binornd(1,p_values(xyz)); % do CCA again with probability p
                CCA_fail(i) = CCA_fail(i)+1; % increment the count of CCA fails
            end
        end
        % X
        % packet_cycle

        if t < K*L % the packet being transmitted may be actual or virtual packet; have to check the packet cycle count
            for n = 1:N % for all nodes
                i = sorted_id(n);
                % i and j are nodes
                for j = [1:i-1 i+1:N] % prev node and next node
                    if tau(j,i) >= 0 %
                        sensed(i,j) = (mod(X(tau(j,i)+1,j),K)==1)*(packet_cycle(j) > 0);
                    else
                        sensed(i,j) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0)*(packet_cycle(j) > 0);
                    end
                end
                sense_flag(i) = sum(sensed(i,1:N)) > 0;
            end
            %        sensed
        else
            for n = 1:N
                i = sorted_id(n);
                for j = [1:i-1 i+1:N]
                    if tau(j,i) >= 0
                        sensed(i,j) = mod(X(tau(j,i)+1,j),K)==1;
                    else
                        sensed(i,j) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0);
                    end
                end
                sense_flag(i) = sum(sensed(i,1:N)) > 0;
            end
        end
        %    sense_flag;

        if t < K*L
            for n = 1:N
                i = sorted_id(n);
                if (mod(X(1,i),K) == 1)&&(t > coll_time(i)-1)
                    for j = [1:i-1 i+1:N]
                        if tau(j,i) >= 0
                            coll_ind(j,i) = ((mod(X(tau(j,i)+1,j),K) == 1) + (mod(X(tau(j,i)+2,j),K) == 1))*(packet_cycle(j) > 0);
                        else
                            coll_ind(j,i) = (X(1,j) == -1)*(sense_flag(j) == 0) ...
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
                            coll_ind(j,i) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0) + (mod(X(1,j),K) == 1);
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
    CCA_fail;
    packet_cycle;
    coll;
    % mean_CCA_fail = mean(CCA_fail);
    % mean_packet_cycles = mean(packet_cycle);
    % mean_coll = mean(coll);

    theta = (sum(packet_cycle)-sum(coll))*K*L/T;
    the(end+1)=theta;
    %print(p_values(j));

    %plot(p_values,theta)
    %end
end
%hold off
the;
%subplot(2,3,5);
plot(p_values,the,'m')
%title("L=" + L + " , D=4")
title("L = 32");
legend('D=0','D=1','D=2','D=3','D=4')
hold off
