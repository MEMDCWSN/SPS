function data = sps_load_bla(filename, type)
% function data = sps_load_bla(filename, type)

fid=fopen(filename,'r');  if fid<=0 fprintf(1,'Error opening file %s!\n',filename); data=[]; return; end;

numLists  = fread(fid,1,'int32');      data = cell(numLists,1);
index = fread(fid,numLists,'int32');   elems = fread(fid,sum(index),type);
elemIdx = 1;
for i=1:numLists
    data{i} = elems(elemIdx:elemIdx+index(i)-1)';
    elemIdx=elemIdx+index(i);
end;

fclose(fid);
