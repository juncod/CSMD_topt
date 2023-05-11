function saveFunction2(saveFileName, x,plotCut)
    forSaveX = x(:);
    saveX = fliplr(reshape(double(forSaveX), [128, 128]))';
    saveMatName = strcat(saveFileName, '.mat');
    saveXlsName = strcat(saveFileName, '.xlsx');
    xlswrite(saveXlsName, saveX)
    save(saveMatName, 'saveX', '-v7.3', '-nocompression')

