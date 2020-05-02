function [ave_sl, std_sl, p20, p15, p10] = ...
    sft_export(PathName, FileName, m_str, m_cov, fillbins, umperpix, ...
    blocksize, scanjump, timer, data, superiorang, OI, AI, CMI, percsarc)

% Exports data arrays into .m files that contain X, Y, distance from edge,
% pattern strength, patern direction, and pattern spacing.  Also exports
% the data into an Excel workbook.

disp(['Exporting data to ' PathName FileName '.m'])
disp(['Exporting data to ' PathName FileName ' Excel Workbook'])

% percsarc = (nnz(m_str(m_cov)>.5)) / (nnz(m_cov)) * 100;

savefile = [PathName 'output_' FileName '.m'];
params.umperpix=umperpix;
params.blocksize=blocksize;
params.scanjump=scanjump;
params.timer=timer;



dataexport = {'x', 'y', 'edgedist(um)', 'strength', 'direction',...
    'sarclength', 'lanewidth', 'OI', 'AI', 'CMI'};
dataexport2 = [data(:,1) data(:,2) data(:,3) ...
    data(:,4) data(:,5) data(:,6) data(:,7) data(:,8) data(:,9) data(:,10)];
% scatter(data((data(:,4)>0.5),3),data((data(:,4)>0.5),6))
% paramexport = {'path', 'file', 'umperpix', 'blocksize', 'scanjump', 'time'; ...
%     PathName FileName params.umperpix params.blocksize params.scanjump ...
%     params.timer};


paramexport = {'path' PathName; 'file' FileName; 'umperpix' params.umperpix; ...
    'blocksize' params.blocksize; 'scanjump' params.scanjump; ...
    'time' params.timer};

warning off MATLAB:xlswrite:AddSheet
xlswrite([PathName FileName '.xls'], dataexport, 'Data Export');
xlswrite([PathName FileName '.xls'], dataexport2, 'Data Export', 'A2');
xlswrite([PathName FileName '.xls'], paramexport, 'Parameters');


superiorang = ceil(superiorang);
fillbins3 = [fillbins;fillbins;fillbins];
sumexport = {'% Area Containing Sarcomeres(>.5)' percsarc; ...
    'Average Sarcomere Strength' OI; ...
    'Superior Angle' superiorang;...
    '% Weighted Alignment Within 20 Degrees' 100*sum(fillbins3(160+superiorang:200+superiorang))/sum(fillbins);... 
'% Weighted Alignment Within 15 Degrees' 100*sum(fillbins3(165+superiorang:195+superiorang))/sum(fillbins);...
'% Weighted Alignment Within 10 Degrees' 100*sum(fillbins3(170+superiorang:190+superiorang))/sum(fillbins); ...
'Average Sarcomere Length' mean(data(:,6)) ;...
'Sarcomere Length Std. Dev.' std(data(:,6));...
'Filtered (strength>0.5) Sarc Length' mean(data(data(:,4)>0.5,6));...
'Filtered (strength>0.5) Sarc Length Std Dev' std(data(data(:,4)>0.5,6));...
'ORGANIZATION INDEX: ' OI ;...
'ALIGNMENT INDEX (MRL) (0 to 1): ' AI ;...
'COMBINED MYOFIBRIL INDEX: ' CMI };
xlswrite([PathName FileName '.xls'], sumexport, 'Summary');

save([PathName FileName '.mat'], 'PathName','FileName','params','data','sumexport');


disp(' ');
disp(['Average Sarcomere Strength: ' num2str(OI) ]);
disp(['% Area Containing Sarcomeres(>.5): ' num2str(percsarc) '%']);
disp(['Superior Angle: ' num2str(superiorang) ' degrees']);
disp(['% Weighted Alignment Within 20 Degrees: ' num2str(100*sum(fillbins3(160+superiorang:200+superiorang))/sum(fillbins)) '%']); 
disp(['% Weighted Alignment Within 15 Degrees: ' num2str(100*sum(fillbins3(165+superiorang:195+superiorang))/sum(fillbins)) '%']);
disp(['% Weighted Alignment Within 10 Degrees: ' num2str(100*sum(fillbins3(170+superiorang:190+superiorang))/sum(fillbins)) '%']);
disp(['Average Sarcomere Length: ' num2str(mean(data(:,6))) ' um' ]);
disp(['Sarcomere Length Std. Dev.: ' num2str(std(data(:,6))) ' um']);
disp(['Average Sarc Length in Strong Pattern Areas (>.5): ' num2str(mean(data(data(:,4)>0.5,6))) ' um']);
disp(['Sarc Length Std. Dev. in Strong Pattern Areas (>.5): ' num2str(std(data(data(:,4)>0.5,6))) ' um']);
disp(' ')
disp('-----------------------------------------')
disp(['            ORGANIZATION INDEX: ' num2str(OI) ]);
disp(['ALIGNMENT INDEX (MRL) (0 to 1): ' num2str(AI) ]);
disp(['      COMBINED MYOFIBRIL INDEX: ' num2str(CMI) ]);
disp('-----------------------------------------')
disp(' ')


percsarc = percsarc/100;
p20 = sum(fillbins3(160+superiorang:200+superiorang))/sum(fillbins);
p15 = sum(fillbins3(165+superiorang:195+superiorang))/sum(fillbins);
p10 = sum(fillbins3(170+superiorang:190+superiorang))/sum(fillbins);
ave_sl = mean(data(:,6));
std_sl = std(data(:,6));
end