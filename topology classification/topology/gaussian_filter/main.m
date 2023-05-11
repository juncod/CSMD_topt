partial_1 = readtable("partial_1.csv");
partial_1_1 = table2array(partial_1);
partial_1_2 = partial_1_1(2:end,2:end);
partial_1_3 = imgaussfilt(partial_1_2,2);

partial_2 = readtable("partial_2.csv");
partial_2_1 = table2array(partial_2);
partial_2_2 = partial_2_1(2:end,2:end);
partial_2_3 = imgaussfilt(partial_2_2,2);

partial_3 = readtable("partial_3.csv");
partial_3_1 = table2array(partial_3);
partial_3_2 = partial_3_1(2:end,2:end);
partial_3_3 = imgaussfilt(partial_3_2,2);

partial_4 = readtable("partial_4.csv");
partial_4_1 = table2array(partial_4);
partial_4_2 = partial_4_1(2:end,2:end);
partial_4_3 = imgaussfilt(partial_4_2,2);


plotCut = 0.2;
savefunction1('partial_gauss_1', partial_1_3, plotCut)
savefunction1('partial_gauss_2', partial_2_3, plotCut)
savefunction1('partial_gauss_3', partial_3_3, plotCut)
savefunction1('partial_gauss_4', partial_4_3, plotCut)

% savefunction2('heat_nocut', heat3, plotCut)
% savefunction2('compl_nocut', comp3, plotCut)