function sps_report(specs, components, minTagLen, specsOri, abinfo, peakTol, reportFilesPath)
% function sps_report(specs, components, minTagLen, specsOri, abinfo, peakTol, reportFilesPath)

if nargin<5 specsOri=[]; end;
if nargin<6 abinfo=[]; end;
if nargin<7 peakTol=0; end;
idxPep = size(specsOri,2);

reportNewFN = 'sps_report.html';
fidNew = fopen(reportNewFN,'w'); if fidNew<=0 fprintf(1,'ERROR opening %s!\n',reportNewFN); end;

if size(specs,2)==5 idxSpec=2; else idxSpec=3; end;
numSpecs = size(specs,1);   peptides = cell(numSpecs,2);   info = zeros(numSpecs,2);   dbInfo = cell(numSpecs,2);
for i=1:numSpecs
    peptides{i,1} = aux_peptide(specs{i,idxSpec},1);   
    info(i,:) = [aux_findLongestTag(peptides{i,1},.5) length(components{i})];
    peptides{i,2} = aux_peptide(specs{i,idxSpec},2);
end
[foo,idxS] = sort(info(:,1));   idxS = idxS(numSpecs:-1:1);   info = [idxS info(idxS,:)];   peptides = peptides(idxS,:);   dbInfo = dbInfo(idxS,:);

idxOk = find(info(:,2)>=minTagLen);
peptides = peptides(idxOk,:);   info = info(idxOk,:);   dbInfo = dbInfo(idxOk,:);

fprintf(fidNew,'<HTML><HEAD><TITLE>Assembly report: %s</TITLE></HEAD><BODY>\n',reportNewFN);
fprintf(fidNew,'<TABLE><TH>Index<TH>Tag Len<TH>Num Specs<TH>Peptide<TH>Peptide expanded (note L = I/L)\n');
numEntries = size(peptides,1);   reportNew=[];
for i=1:numEntries
    fprintf(fidNew,'<TR><TD>%d<TD>%d<TD>%d<TD>%s<TD>%s\n',info(i,1),info(i,2),info(i,3),peptides{i,1},peptides{i,2});
    reportNew = [reportNew; { info(i,1),info(i,2),info(i,3),peptides{i,1},peptides{i,2} }];
end;
fprintf(fidNew,'</TABLE></BODY></HTML>\n');   fclose(fidNew);

function strDenovo = aux_peptide(spec, szJumps, in_masses, in_letters, peakTol)
% function strDenovo = aux_peptide(spec, szJumps, in_masses, in_letters, peakTol)

global AAmasses AAletters;
if nargin<2 szJumps=1; end;
if nargin<4 szJumps=1; in_masses=AAmasses; in_letters=AAletters; end;
if nargin<5 szJumps=1; in_masses=AAmasses; in_letters=AAletters; peakTol=0.5; end;

if szJumps==1 [jumps, letters] = aux_getjumps(in_masses, in_letters, szJumps);
else [jumps, letters] = aux_getjumpsAll(in_masses, in_letters, 2, .01); end;

strDenovo = '';
numPeaks = size(spec,1);   if numPeaks==0 return; end;
if spec(1,1)<=0.0001 masses = spec(2:numPeaks,1)-spec(1:numPeaks-1,1);
else masses = spec(:,1)-[0; spec(1:numPeaks-1,1)]; end;
for m=1:length(masses)
    idx = find(abs(jumps-masses(m))<=peakTol);
    if ~isempty(idx)
        if szJumps==1 aaStr=sprintf('%s',letters{idx}); else aaStr=aux_mergeStrings(letters(idx)); end;
        strDenovo = sprintf('%s%s',strDenovo,aaStr);
    else strDenovo = sprintf('%s[%.1f]',strDenovo,masses(m)); end;
end;

function s = aux_mergeStrings(strs)
% function s = aux_mergeStrings(strs)

single = {};   for i=1:size(strs,1) single = [single; strs{i}]; end;
sizes = zeros(size(single,1),1);   for i=1:size(single,1) sizes(i)=size(single{i},2); end;
single = single(find(sizes==min(sizes)));

if size(single,1)==1 s=single{1}; else 
    s=sprintf('[%s',single{1}); 
    for i=2:size(single,1) str=sprintf('%s,%s',s,single{i}); s=str; end; 
    str=sprintf('%s]',s); s=str;
end;

function [masses, letters] = aux_getjumps(inMasses, inLetters, k)
% function [masses, letters] = aux_getjumps(inMasses, inLetters, k)

if k==0 masses=[]; letters=''; return; end;

baseMasses = inMasses;   resolution = .1;
useMasses = round(inMasses/resolution);    masses = unique(useMasses);    letters = cell(size(masses));   
for iter=1:k
	for i=1:size(useMasses,1)
        idx = find(masses==useMasses(i));
        if iter<2
            if isempty(letters{idx,1})  letters{idx,1} = inLetters(i);
            else
                letters{idx,1} = sprintf('[%s,%s]', letters{idx,1}, inLetters(i));
            end
        else
            letters{idx,1} = sprintf('[%.1f]',masses(idx)/10);
        end;
    end;
    if iter<k
        inMasses = unique(reshape( repmat(inMasses, 1, size(baseMasses,1)) + repmat(baseMasses', size(inMasses,1), 1) , size(inMasses,1)*size(baseMasses,1), 1));
        useMasses = setdiff(unique(round(inMasses/resolution)), masses);
        masses = [masses; useMasses];
    end;
end;
masses = masses * resolution;

function [masses, letters] = aux_getjumpsAll(inMasses, inLetters, k, resolution)
% function [masses, letters] = aux_getjumpsAll(inMasses, inLetters, k, resolution)

if nargin<4 resolution=0.1; end;
if k==0 masses=[]; letters=''; return; end;

baseMasses = inMasses;
useMasses = round(inMasses/resolution);    [masses, idxUnique] = unique(useMasses);    
lettersU = cell(size(idxUnique,1),1);      for i=1:size(idxUnique,1) lettersU{i} = inLetters(idxUnique(i)); end;
maxMass = max(masses);    massesArray = zeros(k*maxMass,1);   lettersArray = cell(k*maxMass,1);

[seqs, massesSeqs] = aux_genAllPerms(masses,lettersU,k);   numSeqs = size(seqs,1);
for i=1:numSeqs
    massesArray(massesSeqs(i)) = massesSeqs(i);   
    lettersArray{massesSeqs(i)}=[lettersArray{massesSeqs(i)};seqs(i)];
end
idx = find(massesArray>0);
masses = massesArray(idx)*resolution;
letters = lettersArray(idx);

function [seqs,massesSeqs] = aux_genAllPerms(masses,letters,k)
% function [seqs,massesSeqs] = aux_genAllPerms(masses,letters,k)

if k==1 seqs=letters; massesSeqs=masses; return; end;   
[seqsTmp,massesTmp] = aux_genAllPerms(masses,letters,k-1);   szLetters = size(letters,1);
seqs = repmat(seqsTmp, szLetters, 1);   massesSeqs = repmat(massesTmp, szLetters, 1);
for i=1:szLetters
    baseIdx = (i-1)*szLetters;
    for j=1:szLetters
        seqs{baseIdx+j}=strcat(letters{i},seqs{baseIdx+j});
        massesSeqs(baseIdx+j)=massesSeqs(baseIdx+j)+masses(i);
    end
end
seqs = [seqsTmp; seqs];   massesSeqs = [massesTmp; massesSeqs];

function tagLen = aux_findLongestTag(tag, peakTol)
% function tagLen = aux_findLongestTag(tag, peakTol)

global AAmasses AAletters;
peakTolInt = round(10*peakTol);   tolRange = [-peakTolInt:peakTolInt];   szTolRange = size(tolRange,2);
jumps1 = unique(repmat(round(10*AAmasses),1,szTolRange)+repmat(tolRange,size(AAmasses,1),1));   jumps1ok = zeros(max(jumps1),1);   jumps1ok(jumps1)=1;   clear jumps1;
jumps2 = aux_getjumps(AAmasses,AAletters,2);   jumps2 = unique(repmat(round(10*jumps2),1,szTolRange)+repmat(tolRange,size(jumps2,1),1));   jumps2ok = zeros(max(jumps2),1);   jumps2ok(jumps2)=1;   clear jumps2;

masses = round(10*aux_getmassesExt(tag,'',[],100));   masses = masses(find(masses>0));
szMasses = size(masses,2);   massType = zeros(1,szMasses);
for i=1:szMasses
    if masses(i)<=size(jumps1ok,1) & jumps1ok(masses(i))==1 massType(i)=1; continue; end;
    if masses(i)<=size(jumps2ok,1) & jumps2ok(masses(i))==1 massType(i)=2; end;
end;

tagLen=0;
tagStart = min(find(massType==1));  if isempty(tagStart) if ~isempty(find(massType==2)) tagLen=2; end; return; end;
tagEnd = tagStart;   jumpPos = 0;
idxJumps = max(find(find(massType==2)<tagStart));  if ~isempty(idxJumps) tagStart=idxJumps; jumpPos=tagStart; end; 
while tagEnd<=szMasses
    if massType(tagEnd)==2 
        if jumpPos>0 tagLen=max(tagLen,tagEnd-tagStart+1); tagStart = jumpPos+1; end;
        jumpPos=tagEnd;
    end
    if massType(tagEnd)==0
        if jumpPos>0 tagLen=max(tagLen,tagEnd-tagStart+1); else tagLen=max(tagLen,tagEnd-tagStart); end;
        tagStart = tagEnd+1;   jumpPos=0;
    end
    tagEnd=tagEnd+1;
end
if jumpPos>0 tagLen=max(tagLen,tagEnd-tagStart+1); else tagLen=max(tagLen,tagEnd-tagStart); end;

function uPeps = aux_uniquePeptides(peptides,mods,modMasses,addModMasses)
% function uPeps = aux_uniquePeptides(peptides,mods,modMasses,addModMasses)

uPeps = unique(peptides);  szU = size(uPeps,1);   uPeps = [uPeps cell(szU,3)];
for i=1:szU
    uPeps{i,3} = find(strcmp(uPeps{i,1},peptides))';
    uPeps{i,2} = size(uPeps{i,3},2);
    if addModMasses>0 | ~isempty(mods) uPeps{i,4} = size(aux_getmassesExt(uPeps{i,1},mods,modMasses,addModMasses),2);
    else uPeps{i,4} = -1; end;
end

function masses = aux_getmasses(peptide);
% function masses = aux_getmasses(peptide);

global AAletters AAmasses;
if isempty(peptide) masses=[]; return; end;
szPeptide = size(peptide,2);
masses = zeros(1,szPeptide);

for i=1:size(AAletters,1)
    idx = find(peptide == AAletters(i));
    masses(idx) = AAmasses(i);
end;


function masses = aux_getmassesExt(peptide,mods,modMasses,addModMasses)
% function masses = aux_getmassesExt(peptide,mods,modMasses,addModMasses)

masses = [];
while ~isempty(peptide)
    if peptide(1)~='['
        newMass = aux_getmasses(peptide(1));
    else
        [strMass, peptide] = strtok(peptide(2:size(peptide,2)),']');
        newMass = str2num(strMass);  % multi-aminoacid jump
        if isempty(newMass) 
            newMass = aux_getmasses(strMass(1));   % equivalent amino acids [I,L] or [Q,K]
        else
            if abs(newMass)<=addModMasses & ~isempty(masses)
                masses(size(masses,2)) = masses(size(masses,2))+newMass;
                peptide = peptide(2:size(peptide,2));
                continue; 
            end;
        end;  
    end;
    peptide = peptide(2:size(peptide,2));
    if ~isempty(peptide)
        p=findstr(peptide(1),mods); if ~isempty(p) newMass = newMass+modMasses(p); peptide = peptide(2:size(peptide,2)); end;
    end
    masses = [masses newMass];
end
