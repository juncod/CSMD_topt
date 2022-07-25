function saveFunction(saveFileName,x,f1,f3,plotCut)
    forSaveX = x(:)>plotCut;
    saveX=fliplr(reshape(double(forSaveX),[100,100]))';
    saveMatName = strcat(saveFileName,'.mat');
    saveXlsName = strcat(saveFileName,'.xlsx');
    saveDensName = strcat(saveFileName,'_dens');
    xlswrite(saveXlsName,saveX)
    save(saveMatName,'saveX','-v7.3','-nocompression')
    saveas(f1,saveDensName,'jpg')
    saveas(f3,saveFileName,'jpg')