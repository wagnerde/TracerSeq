function header = import_csv_header(filename)
% Usage: header = import_csv_header(filename)
% 
% Imports header column from a transcriptome counts csv file.  
%
%
%% CODE:

fid = fopen(filename);
header = cellstr(strsplit(string(textscan(fid,'%s',1)),','))';
header(1)=[];
fclose(fid);

