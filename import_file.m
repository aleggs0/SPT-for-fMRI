function [XX,XXorig] = import_file(subjectno)
%IMPORT_FILE imports rest_mcf file from drive and returns preprocessed and
%original files
%   Detailed explanation goes here
text = fileread('/home/jada2/Data_Store/Warwick_Data/aston1/Brain_Imaging/Beijing_Zang/Beijing_Zang_subjects.txt');
subjectnames = split(text); subjectnames=subjectnames(1:end-1);
destfolder="/mhome/damtp/r/ay343/code/";
destname=destfolder+"rest_mcf.nii";
if ismember(subjectno,[25 30 32 57 107 119])
    disp("warning: the scan for this subject was found to be anomalous; the data may be corrupted")
end
% 128 missing
% 25 skull noise
% 30,32,57 corrupted
% 31 somewhat corrupted
% 119 weird shake
subjectname=string(subjectnames(subjectno));
disp(subjectname)
sourcename="/home/jada2/Data_Store/Warwick_Data/aston1/Brain_Imaging/Beijing_Zang/"+subjectname+"/func/rest_mcf.nii";
if isfile(sourcename)
    copyfile(sourcename,destname)
elseif subjectname=="sub62966"
    return
else
    gunzip(sourcename+".gz",destfolder)
end
[XX, XXorig] = preprocess(destname);
end