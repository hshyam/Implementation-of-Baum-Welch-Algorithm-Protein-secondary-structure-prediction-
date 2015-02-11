function[seqs] = readSeqs(fileName)

% Open the file (train)
FID = fopen(fileName);
line = fscanf(FID,'%s');

% Now find the start and end points of the sequences
p = strfind(line, '<>');
q = strfind(line, 'end');
startEnd = horzcat(p,q);
startEnd = sort(startEnd);
seqs = cell(1,1);
q=0;
for a = 1: length(startEnd)-1
    if (startEnd(a+1)-startEnd(a)>4)
        q=q+1;
        seqs{q,1} = line(startEnd(a)+2:startEnd(a+1)+1);
        seqs{q,1} = strrep(seqs{q,1},'en','');
        seqs{q,1} = strrep(seqs{q,1},'<>','');   
    end
end
    seqs{q,1} = strrep(seqs{q,1},'<','');
end
