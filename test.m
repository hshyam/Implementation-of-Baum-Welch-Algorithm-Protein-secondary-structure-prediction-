function [scores] = test(t, e)

% Open the file (test)
FID = fopen('protein-secondary-structure.test');
line = fscanf(FID,'%s');
chars = [ 'GAVLIPFYWSTCMNQKRHDE' ];
states = [ 'he_'];

% Now find the start and end points of the sequences
p = strfind(line, '<>');
q = strfind(line, '<end>');
testSeqs = cell(length(p),1);
backward=1;
[x,y] = size(e);

% Put the sequences in lines
for a = 1 : length(p) 
    if p(a)>q(backward)
        backward = backward+1;
    end
    testSeqs{a} = line(p(a)+2:q(backward)-1);
end

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
    
    % Decode
    [~, currentSeq] = ismember(currentSeq,chars);
    [~,forward,backward] = frwdbkwd(currentSeq,t,e);
    pStates = forward.*backward;
    finalStr = '';
    for i = 1 : length(pStates)
        if i == 1 
            finalStr = strcat(finalStr, 'h');
        else
            [~, index] = max(pStates(:,i));

            finalStr = strcat(finalStr, states(index));
        end
    end
    
    % Display the text matched state sequences
    disp(currentStates)
    disp(finalStr)
    L = length(finalStr);
    
    % Compute accuracy of matches
    score = 0;
    for i = 1:L
        if currentStates(i)==finalStr(i)
            score = score +1;
        end
    end
    
    % Store the score in an array and then output or display it
    scores(currentI,1) = score/L;
    disp(score/L) 
end

