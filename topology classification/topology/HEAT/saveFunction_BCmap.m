function saveFunction_BCmap(saveBCname, BC_map)
%     saveMatName = strcat(saveBCName,'.mat');
    saveXlsName = strcat(saveBCname,'.xlsx');
    xlswrite(saveXlsName,BC_map)
%     save(saveMatName,'BC_map','-v7.3','-nocompression')
end