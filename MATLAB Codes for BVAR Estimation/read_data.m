% M = readtable("Reg_Data1993.csv");
M = readtable("Reg_Data1981.csv");
Yvars = {'GDPC1', 'CPITOT','UNRATE','FEDFUNDS'}; % subject to change
% Yvars = {'OUTPUTGAP', 'CPITOT','UNRATE','FEDFUNDS'};
% Yvars = {'GDPC1', 'PCE','GPDIC1','FEDFUNDS'}; 
% Yvars = {'GDPC1', 'CPITOT','UNRATE','FEDFUNDS'};
% Yvars = {'GDPC1', 'PCE','GPDIC1','FEDFUNDS'};
Yraw = table2array(M(:, Yvars));
save ('Yraw.dat','Yraw','-ASCII');
W_T = table2array(M(:, {'AVG2_TOT'}));
% W_T = table2array(M(:, {'EPU'}));
save ('W_T.dat','W_T','-ASCII');
W_S = table2array(M(:, {'AVG2_SCORE'}));
save ('W_S.dat','W_S','-ASCII');
W_others = table2array(M(:, {'EPU'}));
save ('W_others.dat','W_others','-ASCII');
yearlab = table2array(M(:, {'yearlab'}));
save ('yearlab.dat','yearlab','-ASCII');