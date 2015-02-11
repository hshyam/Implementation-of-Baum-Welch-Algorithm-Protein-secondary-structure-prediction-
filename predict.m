% Calculate the emission probability
function [emissions, transitions] = predict(seqs)
chars = ['GAVLIPFYWSTCMNQKRHDE'];
states = [ 'he_'];
emissions = zeros(length(states),length(chars));
    for i = 1 : length(seqs)
        currentSeq = seqs{i,1};
        curState = seqs{i,2};
        for p = 1:length(curState)
            emission = currentSeq(p);
            state = curState(p);
            x = strfind(chars,emission);
            y = strfind(states,state);
            emissions(y,x) = emissions(y,x) + 1;
        end 
    end

    % For emissions, find the sum and percentage of probabilities per row 
    rowTotal = sum(emissions,2);
    emissions(1,1:20) = emissions(1,1:20)/rowTotal(1);
    emissions(2,1:20) = emissions(2,1:20)/rowTotal(2);
    emissions(3,1:20) = emissions(3,1:20)/rowTotal(3);

    % Find the transition probability
    transitions =  zeros(length(states));
    for i = 1 : length(seqs)
        statesInput = seqs{i,2};
        for p = 2:length(statesInput)
            lastState = statesInput(p-1);
            curState = statesInput(p);
            x = strfind(states,lastState);
            y = strfind(states,curState);
            transitions(y,x) = transitions(y,x) + 1;
        end   
    end
    
    % For transitions, find the sum and percentage of probabilities per row
    rowTotal = sum(transitions,2);
    transitions(1,1:y) = transitions(1,1:y)/rowTotal(1);
    transitions(2,1:y) = transitions(2,1:y)/rowTotal(2);
    transitions(3,1:y) = transitions(3,1:y)/rowTotal(3);
 end