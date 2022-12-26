function [NODE,ELEM] = inp_(filename)
fid = fopen(filename);
key = 0;
nodekey = 0;
elemkey = 0;

while key == 0 %% Find Start point of NODE
    tline = fgetl(fid);
    if strcmp(tline, '*Node')
        key = 1;
        break;
    end
end

NODE = [];
while nodekey == 0 %% Import NODE data & Find Start point of ELEM
    tline1 = fgetl(fid);
    if strcmp(tline1, '*Element, type=CPS4R')
        nodekey = 1;
        break;
    else
        node = str2num(tline1);
        NODE = [NODE; node];
    end
end

ELEM = [];
while elemkey == 0 %% Import ELEM data
    tline2 = fgetl(fid);
    if strcmp(tline2, '*End Part')
        elemkey = 1;
        break;
    else
        elem = str2num(tline2);
        ELEM = [ELEM; elem];
    end
end

NODE(:,1) = [];  % Remove first column
ELEM(:,1) = [];

end