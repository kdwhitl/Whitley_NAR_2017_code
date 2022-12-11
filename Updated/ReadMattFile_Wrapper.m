function data = ReadMattFile_Wrapper(startpath,rootfile)
%Call this function first, it will then call the appropriate data reading
%function

%before adding version numbers
if str2double(rootfile(1:6)) < 90507
    data = ReadMattFile_TimeSharing(startpath,rootfile);
else
%     data = ReadMattFile_v8(startpath,rootfile);
    data = ReadJointFile(startpath,rootfile);
end
