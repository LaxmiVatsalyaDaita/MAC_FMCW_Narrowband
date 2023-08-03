%single CCA throughput calc
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
    %CCA2_fail = zeros(1,N);
    coll_ind = zeros(N,N);
    coll =  zeros(1,N);
    packet_cycle = zeros(1,N);
    sense_flag =  zeros(1,N);
    %sense2_flag = zeros(1,N);
    sensed = zeros(N);
    %sensed2 = zeros(N);
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

    coll_time = K*L-X(1,:);
    for t = 1:T
        X(2:D+2,:) = X(1:D+1,:);
        for i = 1:N
            if (X(2,i) == K*L)||(X(2,i) == 0)
                X(1,i) = -1*binornd(1,p_values(xyz));
            end
            if (X(2,i) > 0)&&(X(2,i) < K*L)
                X(1,i) = X(2,i)+1;
            end
            if (X(2,i) == -1) && (sense_flag(i) == 0)
                X(1,i) = 1;
                packet_cycle(i) = packet_cycle(i) + 1;
            end
            if (X(2,i) == -1) && (sense_flag(i) == 1)
                X(1,i) = -1*binornd(1,p_values(xyz));
                CCA_fail(i) = CCA_fail(i)+1;
            end
        end
        % X
        % packet_cycle

        if t < K*L
            for n = 1:N
                i = sorted_id(n);
                for j = [1:i-1 i+1:N]
                    if tau(j,i) >= 0
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
                            coll_ind(j,i) = ((mod(X(tau(j,i)+1,j),K) == 1))*(packet_cycle(j) > 0);
                        else
                            coll_ind(j,i) = (X(1,j) == -1)*(sense_flag(j) == 0) ...
                                + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0)*(packet_cycle(j) > 0);
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
                            coll_ind(j,i) = (mod(X(tau(j,i)+1,j),K) == 1);
                        else
                            coll_ind(j,i) = (X(1,j) == -1)*(sense_flag(j) == 0) + (X(1,j) > 0)*(X(1,j) < K*L)*(mod(X(1,j),K) == 0);
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
    theta = (sum(packet_cycle)-sum(coll))*K*L/(T*40);
    the(end+1)=theta;
end

the