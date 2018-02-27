function [vSets, eSets] = sps_assembleE(aligns)
% function [vSets, eSets] = sps_assembleE(aligns)

numPeptides = size(unique(aligns(:,[1 2])),1);
vSets = cell(numPeptides,2);
contigsIdx = zeros(max(max(aligns(:,[1 2]))),1);
numContigSets = 0;

szAligns = size(aligns,1);
toProcess = ones(szAligns,1); idxAligns = [1:szAligns]';   tensLeft = floor(szAligns/10000);
for idx=1:szAligns
    
    i = aligns(idx,1);
    j = aligns(idx,2);
    
    if contigsIdx(i)>0 & contigsIdx(j)>0 & contigsIdx(i)==contigsIdx(j) 
        continue; 
    end;  
    
    idxI = contigsIdx(i);
    idxJ = contigsIdx(j);
    
    if idxI==0 & idxJ==0 numContigSets=numContigSets+1; vSets{numContigSets,1} = [i j];   contigsIdx([i j])=numContigSets;   
    else if idxI==0 vSets{idxJ,1} = [vSets{idxJ,1} i];  contigsIdx(i)=idxJ;  
        else if idxJ==0 vSets{idxI,1} = [vSets{idxI,1} j];  contigsIdx(j)=idxI;
            else if idxI~=idxJ 
                    vSets{idxI,1} = unique([vSets{idxI,1} vSets{idxJ,1}]); 
                    contigsIdx(vSets{idxJ,1}) = idxI;
                    if idxJ<numContigSets
                        vSets{idxJ,1} = vSets{numContigSets,1}; 
                        contigsIdx(vSets{numContigSets,1}) = idxJ;
                    end
                    numContigSets=numContigSets-1;
                end;
            end;
        end;
    end;

    if max(contigsIdx)>numContigSets 
        a=1;
    end
    
    curTensLeft=floor(size(idxAligns,1)/10000); 
end; % while
vSets = vSets(1:numContigSets,:);
edgesSet = contigsIdx(aligns(:,1));   eSets = cell(size(vSets,1),1);
for i=1:size(vSets,1)  eSets{i} = find(edgesSet==i)'; end;

function idx = findContig(vertex, vSets)

idx=0;
for c=1:size(vSets,1)
    if ismember(vertex, vSets{c,1}) idx=c; return; end;
end;
