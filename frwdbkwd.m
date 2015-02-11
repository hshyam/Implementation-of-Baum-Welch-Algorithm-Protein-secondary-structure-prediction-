function[probabilitySeq, forward, backward, probability] = frwdbkwd(seq,tr,e)

% Variable initialization
len = length(seq);
states = size(tr,1);
forward = zeros(states,len);
forward(1,1) = 1;
probability = zeros(1,len);
probability(1) = 1;

    % Find forward probabilities 
    for i = 2:len
        for state = 1:states
        sumC = sum(forward(:,i-1) .*tr(:,state));
        forward(state,i) = e(state,seq(i)) .* sumC;
    end
    
    % Prevent overflow by summing each column to 1 
    probability(i) = norm(forward(:,i),1);
    forward(:,i) =  forward(:,i)./probability(i);
    end

    % Find backward probabilities 
    backward = ones(states,len);
        for i = len-1:-1:1
            for state = 1:states
                 temp = tr(state,:)'.* backward(:,i+1) .* e(:,seq(i+1));
                 sumColumn = sum(temp);
                 backward(state,i) = sumColumn/probability(i+1) ; 
            end
        end

% Computation of the sequence probability
probabilitySeq = sum(log(probability));




