N=50;
K=40;
T=10000;
D=2;
L=1;
p_values = [0:0.1:1];

X_(1,:) = randi([1 K*L],1,N);
del_ = rand(1,N); %clock offsets with respect to an universal clock

tri = triu(rand(N));
the = [];
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
the
%figure;
%subplot(2,3,1);
plot(p_values,the,'r')
%title("L=" + L + " , D=0")
hold on


% 2 CCA and L = 1

the = [];
%function theta = CSMA_noncollocated_twocca(N,K,L,p,T)

for xyz = 1:length(p_values)

    CCA1_fail = zeros(1,N);
    CCA2_fail = zeros(1,N);
    coll_ind = zeros(N,N);
    coll =  zeros(1,N);
    packet_cycle = zeros(1,N);
    sense1_flag =  zeros(1,N);
    sense2_flag = zeros(1,N);
    sensed1 = zeros(N);
    sensed2 = zeros(N);
    %D = 0;
    X = zeros(D+2,N);
    %X(1,:) = randi([1 K*L],1,N);
    %del = rand(1,N); %clock offsets with respect to an universal clock
    X(1,:) = X_(1,:);
    del = del_;
    [~, sorted_id] = sort(del);
    delta = repmat(del',1,N) - repmat(del,N,1); %relative clock offset matrix
    %u = D*triu(rand(N));
    u = D*tri;
    d = u + u'- 2*diag(u(1:N+1:N^2));
    tau = floor(d + delta);
    %Random states correspond to initial dummy packet transmissions.
    %Actual activities start at random times, at the end of dummy packets.

    coll_time = K*L-X(1,:); %This is to avoid overcounting of collisions
    for t = 1:T
        X(2:D+2,:) = X(1:D+1,:);
        for i = 1:N
            if (X(2,i) == K*L)||(X(2,i) == 0)
                X(1,i) = -1*binornd(1,p_values(xyz));
            end
            if (X(2,i) > 0)&&(X(2,i) < K*L)
                X(1,i) = X(2,i)+1;
            end
            if (X(2,i) == -1) && (sense1_flag(i) == 0)
                X(1,i) = -2;
            end
            if (X(2,i) == -1) && (sense1_flag(i) == 1)
                X(1,i) = -1*binornd(1,p_values(xyz));
                CCA1_fail(i) = CCA1_fail(i)+1;
            end
            if (X(2,i) == -2) && (sense2_flag(i) == 0)
                X(1,i) = 1;
                packet_cycle(i) = packet_cycle(i) + 1;
            end
            if (X(2,i) == -2) && (sense2_flag(i) == 1)
                X(1,i) = -1*binornd(1,p_values(xyz));
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
                    %    sensed1(i,j) = (X(tau(j,i)+1,j) > 0)*(X(tau(j,i)+1,j) <= (K-1)*L)*((mod(X(tau(j,i)+1,j),K)==K-1)+(mod(X(tau(j,i)+1,j),K) == 0))*(packet_cycle(j) > 0);
                        sensed1(i,j) = (X(tau(j,i)+1,j) > 0)*(mod(X(tau(j,i)+1,j),K)==K-1)*(packet_cycle(j) > 0);
                   
                    else
                      %  sensed1(i,j) = (X(1,j) > 0)*(X(1,j) <= (K-1)*L)*((mod(X(1,j),K)==K-2)+(mod(X(1,j),K) == K-1))*(packet_cycle(j) > 0);
                        sensed1(i,j) = (X(1,j) > 0)*(mod(X(1,j),K)==K-2)*(packet_cycle(j) > 0);   
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
                        %sensed1(i,j) = (X(tau(j,i)+1,j) > 0)*(X(tau(j,i)+1,j) <= (K-1)*L)*((mod(X(tau(j,i)+1,j),K)==K-1)+(mod(X(tau(j,i)+1,j),K) == 0));
                        sensed1(i,j) = (X(tau(j,i)+1,j) > 0)*(mod(X(tau(j,i)+1,j),K)==K-1);
                    else
                        sensed1(i,j) = (X(1,j) > 0)*(mod(X(1,j),K)==K-2);
                        %sensed1(i,j) = (X(1,j) > 0)*(X(1,j) <= (K-1)*L)*(mod(X(1,j),K)==K-2)+(mod(X(1,j),K) == K-1);
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
    theta = (sum(packet_cycle)-sum(coll))*K*L/T;
    the(end+1)=theta;
end
the
%figure;
%subplot(2,3,1);
plot(p_values,the,'c')
%title("L=" + L + " , D=0")




% L = 32

N=50;
K=40;
T=10000;
D=2;
L=28;
p_values = [0:0.1:1];

X_(1,:) = randi([1 K*L],1,N);
del_ = rand(1,N); %clock offsets with respect to an universal clock

tri = triu(rand(N));
the = [];
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
the
%figure;
%subplot(2,3,1);
plot(p_values,the,'y')
%title("L=" + L + " , D=0")



% 2 CCA and L = 32

the = [];
%function theta = CSMA_noncollocated_twocca(N,K,L,p,T)

for xyz = 1:length(p_values)

    CCA1_fail = zeros(1,N);
    CCA2_fail = zeros(1,N);
    coll_ind = zeros(N,N);
    coll =  zeros(1,N);
    packet_cycle = zeros(1,N);
    sense1_flag =  zeros(1,N);
    sense2_flag = zeros(1,N);
    sensed1 = zeros(N);
    sensed2 = zeros(N);
    %D = 0;
    X = zeros(D+2,N);
    %X(1,:) = randi([1 K*L],1,N);
    %del = rand(1,N); %clock offsets with respect to an universal clock
    X(1,:) = X_(1,:);
    del = del_;
    [~, sorted_id] = sort(del);
    delta = repmat(del',1,N) - repmat(del,N,1); %relative clock offset matrix
    %u = D*triu(rand(N));
    u = D*tri;
    d = u + u'- 2*diag(u(1:N+1:N^2));
    tau = floor(d + delta);
    %Random states correspond to initial dummy packet transmissions.
    %Actual activities start at random times, at the end of dummy packets.%del = [0.9842    0.6284    0.5694    0.8980    0.2170    0.4578    0.6538    0.1917    0.7327    0.9100]

    coll_time = K*L-X(1,:); %This is to avoid overcounting of collisions
    for t = 1:T
        X(2:D+2,:) = X(1:D+1,:);
        for i = 1:N
            if (X(2,i) == K*L)||(X(2,i) == 0)
                X(1,i) = -1*binornd(1,p_values(xyz));
            end
            if (X(2,i) > 0)&&(X(2,i) < K*L)
                X(1,i) = X(2,i)+1;
            end
            if (X(2,i) == -1) && (sense1_flag(i) == 0)
                X(1,i) = -2;
            end
            if (X(2,i) == -1) && (sense1_flag(i) == 1)
                X(1,i) = -1*binornd(1,p_values(xyz));
                CCA1_fail(i) = CCA1_fail(i)+1;
            end
            if (X(2,i) == -2) && (sense2_flag(i) == 0)
                X(1,i) = 1;
                packet_cycle(i) = packet_cycle(i) + 1;
            end
            if (X(2,i) == -2) && (sense2_flag(i) == 1)
                X(1,i) = -1*binornd(1,p_values(xyz));
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
                        sensed1(i,j) = (X(tau(j,i)+1,j) > 0)*(mod(X(tau(j,i)+1,j),K)==K-1)*(packet_cycle(j) > 0);
                        %sensed1(i,j) = (X(tau(j,i)+1,j) > 0)*(X(tau(j,i)+1,j) <= (K-1)*L)*((mod(X(tau(j,i)+1,j),K)==K-1)+(mod(X(tau(j,i)+1,j),K) == 0))*(packet_cycle(j) > 0);
                    else
                        sensed1(i,j) = (X(1,j) > 0)*(mod(X(1,j),K)==K-2)*(packet_cycle(j) > 0);
                        %sensed1(i,j) = (X(1,j) > 0)*(X(1,j) <= (K-1)*L)*((mod(X(1,j),K)==K-2)+(mod(X(1,j),K) == K-1))*(packet_cycle(j) > 0);
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
                        sensed1(i,j) = (X(tau(j,i)+1,j) > 0)*(mod(X(tau(j,i)+1,j),K)==K-1);
                        %sensed1(i,j) = (X(tau(j,i)+1,j) > 0)*(X(tau(j,i)+1,j) <= (K-1)*L)*((mod(X(tau(j,i)+1,j),K)==K-1)+(mod(X(tau(j,i)+1,j),K) == 0));
                    else
                        sensed1(i,j) = (X(1,j) > 0)*(mod(X(1,j),K)==K-2);
                        %sensed1(i,j) = (X(1,j) > 0)*(X(1,j) <= (K-1)*L)*(mod(X(1,j),K)==K-2)+(mod(X(1,j),K) == K-1);
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
    theta = (sum(packet_cycle)-sum(coll))*K*L/T;
    the(end+1)=theta;
end
the
%figure;
%subplot(2,3,1);
plot(p_values,the,'m')
%title("L=" + L + " , D=0")
ylim([0 40])


legend('L=1, 1CCA', 'L=1, 2CCA', 'L=28, 1CCA', 'L=28, 2CCA')
title('D=2')

hold off