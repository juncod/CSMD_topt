function saveFunction(saveFileName,x,f3,plotCut)
    forSaveX = x(:)>plotCut;
    saveX=fliplr(reshape(double(forSaveX),[sqrt(length(x)),sqrt(length(x))]))';
    saveMatName = strcat(saveFileName,'.mat');
    saveXlsName = strcat(saveFileName,'.xlsx');
    % saveDensName = strcat(saveFileName,'_dens');
    xlswrite(saveXlsName,saveX)
    save(saveMatName,'saveX','-v7.3','-nocompression')
    % saveas(f1,saveDensName,'jpg')
    saveas(f3,saveFileName,'jpg')