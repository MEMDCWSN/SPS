function abinfo = sps_load_abinfo(filename)
% function abinfo = sps_load_abinfo(filename)

fid = fopen(filename,'r'); if fid<=0 fprintf(1,'Error opening %s!\n',filename); abinfo={}; return; end;

numUsedSpectra = fread(fid,1,'int32');   
if numUsedSpectra>0  % File version 1.0
    specIndices = fread(fid,numUsedSpectra,'int32');
    specInfo = fread(fid,numUsedSpectra*2,'int16');    specInfo = reshape(specInfo, 2, numUsedSpectra)';   
    specInfo = [specIndices+1 specInfo(:,1)+1 specInfo(:,2)];  % +1 to convert all indices to 1-based
    clear specIndices;
else
    vMajor = fread(fid,1,'ushort');   vMinor = fread(fid,1,'ushort');
    if vMajor>1 | (vMajor==1 & vMinor>1)
        fprintf(1,'ERROR (sps_load_abinfo): Cannot read file format version %d.%d\n',vMajor,vMinor);
        return;
    end
    numUsedSpectra = fread(fid,1,'int32');
    specIndices = fread(fid,2*numUsedSpectra,'int32') + 1;   % +1 to convert all indices to 1-based
    specInfo = fread(fid,numUsedSpectra,'int8');
    specInfo = [reshape(specIndices, 2, numUsedSpectra)' specInfo];
    clear specIndices;
end

numComponents = fread(fid,1,'int32');
vertsPerComp = fread(fid,numComponents,'int16');   totNumVerts = sum(vertsPerComp);
peaksPerVert = fread(fid,totNumVerts,'int16');     totNumPeaks = sum(peaksPerVert);
vertsSpecIdx = fread(fid,totNumPeaks,'int32') + 1;   % +1 to convert all indices to 1-based
vertsPeakMass = fread(fid,totNumPeaks,'float32');
fclose(fid);

abinfo = cell(numComponents,2);   curV = 1;   curP = 1;
for c=1:numComponents
    idx = find(specInfo(:,2)==c);
    abinfo{c,1} = specInfo(idx,[1 3]);
    abinfo{c,2} = cell(vertsPerComp(c),1);
    for v=1:vertsPerComp(c)
        peakRange = [curP:curP+peaksPerVert(curV)-1];    curP = curP+peaksPerVert(curV);   curV = curV+1;
        abinfo{c,2}{v} = [vertsSpecIdx(peakRange) vertsPeakMass(peakRange)];
    end
end
