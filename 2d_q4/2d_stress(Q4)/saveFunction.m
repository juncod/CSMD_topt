function saveFunction(saveFileName,x,f1,f3,plotCut)
    forSaveX = x(:)>plotCut;
    saveX=fliplr(reshape(forSaveX,[100,100]))';
    saveMatName = strcat(saveFileName,'.mat');
    saveXlsName = strcat(saveFileName,'.xls');
    saveDensName = strcat(saveFileName,'_dens');
    save(saveXlsName,'saveX')
    save(saveMatName,'saveX','-v7.3','-nocompression')
    saveas(f1,saveDensName,'jpg')
    saveas(f3,saveFileName,'jpg')