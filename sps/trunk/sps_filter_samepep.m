function [idxKeep, info] = sps_filter_samepep(parentMasses, aligns, tolerance)
% function [idxKeep, info] = sps_filter_samepep(parentMasses, aligns, tolerance)

[vSets, eSets] = sps_assembleE(aligns);

numSets = size(vSets,1);   keepSet = zeros(numSets,1);
if size(parentMasses,2)>size(parentMasses,1) parentMasses=parentMasses'; end;
for setIdx=1:numSets
    setShifts = aligns(eSets{setIdx},3:4);
    idxOk = find(abs(setShifts)>tolerance);
    if ~isempty(idxOk) keepSet(setIdx)=1; continue; end;  
    
    setPMs = parentMasses(vSets{setIdx});
    idxOk = find(abs( repmat(setPMs,1,length(setPMs)) - repmat(setPMs',length(setPMs),1) )>tolerance);
    if ~isempty(idxOk) keepSet(setIdx)=1; end;
end

idxKeep = [eSets{find(keepSet==1)}]';
info = [numSets size(find(keepSet==1),1) size(idxKeep,1)];
