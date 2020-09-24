% saves X, T, colors, in the densest naive format I can think of. 
% X is half-ints that are 2 bytes each
% T is singles that are 4 bytes each
% colors is uint8s that are 1 bytes each
% result is written in binary to file so there's no overhead. 
% this is probably not an ideal format to load from.
function SaveConstantTriangulationBinary(outpicname, X, T, colors)
    assert(numel(T)==numel(colors));
    
    T = uint16(T);
    colors = uint8(colors);
    X = single(X);
    nX = uint16(size(X,1));
    nT = uint16(size(T,1));
    
    nXt = typecast(nX,'uint8');
    nTt = typecast(nT,'uint8');
    Xt = typecast(X(:),'uint8');
    Tt = typecast(T(:),'uint8');
    colorst = typecast(colors(:),'uint8');
    
    all = [nXt';nTt';Xt;Tt;colorst];
    fid = fopen(outpicname,'w');
    fwrite(fid, all);
    fclose(fid);

%     expectedSize = 2*2 + nX*2*4 + nT*3*2 + nT*3*1;
%     s = dir(outpicname); actualSize = s.bytes;
%     [expectedSize actualSize];
end