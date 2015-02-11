clear
seqs = readSeqs('protein-secondary-structure.train');
seqs = parse(seqs);
oldsequence = seqs
chars = [ 'GAVLIPFYWSTCMNQKRHDE' ];
states = [ 'he_'];

% Find states and emissions
states = size(states,2);
emissions = size(chars,2);

% Obtain Maximum Likelihood Learning/Estimation
[e, t] = predict(seqs);
hmm_ml{1} = e;
hmm_ml{2} = t;

% Create new sequences from characters 
numOfSeqs = length(seqs);
newSeqs = cell(numOfSeqs,1);
for count = 1:numOfSeqs
    [~, newSeqs{count}] = ismember(seqs{count},chars);
end
seqs = newSeqs
totalprobability = 1;
    while true
    prevTotal = totalprobability;
    totalprobability = 1;
    newE = zeros(states,emissions);
    newTr = zeros(size(t));

      % Obtain a sequence to train and then find the log 
      % of Tr and E to prevent any underflows
      for seqN = 1:length(seqs)
        seq = seqs{seqN};
        seqLength = length(seq);  
        logE = log(e);
        logT = log(t);
        
        % Find probabilites - forward and backward
        [probabilitySeq,forward,backward,probability] = frwdbkwd(seq,t,e);
         forward = log(forward);
         backward = log(backward);
         totalprobability = totalprobability + probabilitySeq;
             
         % Now, update each of the transitions based on the
         % sum of new probabilities
         for p = 1:states
                for q = 1:states
                 temp = 0;
                    for r = 1:seqLength-1
                     sumofLogs = forward(p,r) + logT(p,q) + logE(q,seq(r+1)) + backward(q,r+1);
                     temp = temp + (exp(sumofLogs)/probability(r+1));
                    end
                     newTr(p,q) = newTr(p,q) + temp;
                end
         end
         
         % Next, update of the columns in the sequence where
         % the current iter is the emission 
         for p = 1:states
             for r = 1:emissions
                index = find(seq == r);
                newE(p,r) = newE(p,r) + sum(exp(forward(p,index)+backward(p,index)));
             end
          end
       end
    
    % Matrices updation 
    totalE = sum(newE,2);
    for r = 2 : emissions
        totalE(:,r) = totalE(:,r-1);
    end
    
    e = newE./totalE;
    totalTr = sum(newTr,2);
    
    for r = 2 : states
        totalTr(:,r) = totalTr(:,r-1);
    end
    t  = newTr./totalTr;
    
    % Model precision 
    if abs((totalprobability-prevTotal)/prevTotal) < 1*10^-5
        break
    else
        totalprobability
    end 
end

hmm_bw{1} = e;
hmm_bw{2} = t;

% Run Maximum Liklihood and Baum-Welch tests 
scores_ml = test(hmm_ml{2},hmm_ml{1})
scores_bw = test(hmm_bw{2},hmm_bw{1})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Testing - now find the start and end points of the sequences
p = strfind(line, '<>');
q = strfind(line, '<end>');
testSeqs = cell(length(p),1);
backward=1;
[x,y] = size(e);

% Put the states and sequences in their respective lines
for a = 1 : length(p) 
    if p(a)>q(backward)
        backward = backward+1;
    end
    testSeqs{a} = line(p(a)+2:q(backward)-1);
end

% Finding the most likely sequence of hidden states
for currentI = 1: length(testSeqs)
    current = testSeqs{currentI};
    current = strrep(current, '<>', '');
    L=length(current);
    
    % Create a single string of emissions by stripping out the states
    currentSeq = '';
    for i=1:2:L
        currentSeq = strcat(currentSeq,current(i));
    end
    
    % Now, create a single string of states by stripping out the emissions
    currentStates = '';
    for i=2:2:L
        currentStates = strcat(currentStates,current(i));
    end
    
    % Array initialization
    L = length(currentStates);
    if length(currentStates)~=length(currentSeq)
        pause
    end
    V=zeros(x,L); 
    pi=zeros(1,L);

    % First column of arrray initialization 
    for q=1:x
       V(q,1) = log2(1/x) + E(q,index(currentSeq(1),chars));
    end   	
    
    % Probability score computation for all emissions and states
    for i=2:L
        for r=1:x
           for q=1:x
             temp(q) = V(q,i-1) + A(q,r);
           end
    
          [V(r,i), ptr(r,i)] = max(temp);
          
          V(r,i) = V(r,i) + E(r,index(currentSeq(i),chars));
        end
    end
    
    % Now, we initialize back tracing by finding the last column 
    for q=1:x
       temp(q)=V(q,L);
    end
    [maxscore, pi(L)]=max(temp);

    % Find most likely sequence by back tracing 
    for i=L:-1:2
       pi(i-1) = ptr(pi(i),i);
    end
    
    % Display the text matched state sequences
    probpath = states(pi);
    disp(currentStates)
    disp(probpath)
    
    % Compute accuracy of matches 
    score = 0;
    for i = 1:L
        if currentStates(i)==probpath(i)
            score = score +1;
        end
    end
    
    % Store the score in an array and then output or display it
    scores(currentI,1) = score/L;
    disp(score/L)
end