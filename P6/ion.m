function [ION_ALPHA,ION_BETA]=ion(filename)
fid = fopen(filename);
if fid == -1
    errordlg(['The file ''' filename ''' does not exist.']);
    return;
end

ION_ALPHA  = 0;
while ION_ALPHA == 0
    current_line = fgetl(fid);
    if contains(current_line,'ION ALPHA')
        ION_ALPHA=str2num(current_line(1:60));
    end
end

ION_BETA = 0;
while ION_BETA == 0
    current_line = fgetl(fid);
    if contains(current_line,'ION BETA')
        ION_BETA=str2num(current_line(1:60));
    end
end