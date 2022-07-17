function saveFunction(saveFileName,x)
    forSaveX = x(:)>0.1;
    saveX=fliplr(reshape(forSaveX,[100,100]))';
    saveMatName = strcat(saveFileName,'.mat');
    saveXlsName = strcat(saveFileName,'.xls');
    save(saveXlsName,'saveX')
    save (saveMatName,'saveX','-v7.3','-nocompression')