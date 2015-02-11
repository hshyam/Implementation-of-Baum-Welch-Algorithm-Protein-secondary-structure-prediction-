function[seqs] = parse(seqs)

    % Create a single string of emissions by stripping out the states  
    for i = 1 : length(seqs)
    current = seqs{i};
    L=length(current);
    seq = '';
    for p=1:2:L
        seq = strcat(seq,current(p));
    end
    
    % Now, create a single string of states by stripping out the emissions
    states = '';
    for p=2:2:L
        states = strcat(states,current(p));
    end
        seqs{i,1} = seq;
        seqs{i,2} = states;
    end
end