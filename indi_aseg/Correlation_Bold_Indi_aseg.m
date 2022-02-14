%%
imgDir=dir('*.*');
sheet=2;
[ Num, Txt ] = xlsread( 'D:\jingzongfmri\stimTFnew.xlsx',sheet );
filename = Txt(2:46,[3,4]);
subjnum = length(filename);
Visit = {'visit1','visit2'};
%% FC excle
% excle for pre and post FC
%FCfile='\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\tbl\fMRI_individual_FC_aseg.xlsx';
FCfile='\\172.16.4.19\eeg\Projects\tDCS\202201-SMHC_cathode_FP2\tbl\fMRI_individual_FC_aseg_new.xlsx';

%% set session 

ROIs = {'2012','2014','1027','2027','1028','2028','26','58','1008','2008'};
ROIsName ={'rl-OFC','rm-OFC','ll-PFC','rl-PFC','lm-PFC','rm-PFC','l-NAc','r-NAc','l-IPL','r-IPL'};
FCpath1 ='D:\jingzongfmri\data\fmri_';
FCpath2 ='_individual\Results\ROISignals_FunImgARCF\';
%la:lateral 外侧；me:medial 内侧 IPL :inferiorparietal
%%  write FC （post-pre） to excle
for roi1 = 1:length(ROIs)
    for roi2 = roi1+1:length(ROIs)
        Sheetname = ['FC_',char(ROIsName(roi1)),'_',char(ROIsName(roi2))];
        ColNames = {'index','visit1','visit2','D-value','D-ratio','absD-value','absD-ratio'};
        Table  = table();
        FCRoisAllSubj = zeros(subjnum,length(Visit));
        for subj = 1:subjnum
            for visit=1:length(Visit)
                fcpath = [FCpath1,char(Visit(visit)),FCpath2,'ROICorrelation_FisherZ_',char(filename(subj,1)),'.mat'];
                if exist(fcpath,'file')                  
                    ROIFile = [FCpath1,char(Visit(visit)),FCpath2,'ROI_OrderKey_',char(filename(subj,1)),'.tsv'];
                    fid = fopen(ROIFile);
                    C = textscan(fid, '%s %s %s ', 'HeaderLines', 1);
                    fclose(fid);
                    ROIsList = [C{1,1},C{1,2}];
                    idx1 =  str2num(ROIsList{find(strcmp(ROIsList(:,2), ROIs(roi1))),1}); 
                    idx2 =  str2num(ROIsList{find(strcmp(ROIsList(:,2), ROIs(roi2))),1});  
                    load(fcpath);
                    fc = ROICorrelation_FisherZ;
                    FCRoisAllSubj(subj,visit) = fc(idx1,idx2) 
                end
          
            end
        end
        different = FCRoisAllSubj(:,2)-FCRoisAllSubj(:,1);
        absdifferent = abs(FCRoisAllSubj(:,2))-abs(FCRoisAllSubj(:,1));
        differentR = different./FCRoisAllSubj(:,1);
        absdifferentR = absdifferent./abs(FCRoisAllSubj(:,1));
        Table=table([1:subjnum]',FCRoisAllSubj(:,1),FCRoisAllSubj(:,2),different,differentR,absdifferent,absdifferentR,'VariableNames',ColNames);
        writetable(Table,FCfile,'Sheet',Sheetname)
    end
end



%% Multivarables linear regression for BOLD-ef

clear;
clc;
close all;

dn_results = '\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\CorrBold_individual_aseg\'
dn_root = '\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\tbl\';
dn_EF = [dn_root,'EF.xlsx'];
dn_BOLD =  [dn_root,'fMRI_individual_FC_aseg.xlsx'];
dn_scale = [dn_root,'scales.xlsx'];
dn_SubjInfo =  [dn_root,'BOLD_1st.xlsx'];
[Type,Sheet,Format]=xlsfinfo(dn_BOLD);
calcs_x = {'rl-OFC','rm-OFC','rl-PFC','rm-PFC'};
calcs_y= Sheet;
calcs_yy = {'PSQI','YBOCS','YBOCS_Obsessions','YBOCS_Compulsions','BDI','BAI'};

diffratio = {'D-value','D-ratio'};
abs={'','abs'};
stims = {'stim','sham'};
labels_stim = [1 0];

[ Num, Txt ] = xlsread( dn_SubjInfo,'SubjInfo' );
tbl_info =Num(:,5:7);
temp = find(tbl_info(:,2)==1&tbl_info(:,3)==1);
StimSubjs = find(tbl_info(temp,1)==1 );
ShamSubjs = find(tbl_info(temp,1)==0 );
colorline = lines(length(calcs_x));
X_para={'mean_EF','perc95_EF'};
Y1_para=diffratio ;
Y2_para=abs;
Y_para={'D-value',	'D-ratio',	'absD-value',	'absD-ratio'};
YY_para=diffratio ;
%% sperate stim and sham EF-bold
% Determination of ROI
%if ~exist(fullfile(fileparts(dn_results),'stat','MLR.txt'),'file') mkdir(fullfile(fileparts(dn_results),'stat','MLR.txt'));end
diary(fullfile(fileparts(dn_results),'stat','EF_individual_diffFC_asegMLR.txt'));
for x_para =  1:length(X_para)
    tbl_EF = readtable(fullfile(dn_EF),'Sheet',X_para{x_para},'VariableNamingRule','preserve');
for y_para = 1:length(Y_para)
for i_calc_x = 1:length(calcs_x)
for i_fc = 1:length(calcs_y )
    uy = ''; if mod(y_para,2),else uy = '%'; end
    if y_para>2,  uy = ['abs',uy];end
    label_respondx = [ calcs_x{i_calc_x} '-' X_para{x_para} ];
    label_respondx(find(label_respondx=='_'))='-' ;
    Sheet = char(calcs_y{i_fc});
    tbl_outcomes = readtable(fullfile(dn_BOLD),'Sheet',Sheet,'VariableNamingRule','preserve'); 
    label_respondy = [ '\Delta' calcs_y{i_fc} ')' uy ];
    label_respondy(find(label_respondy=='_'))='(&' ;
    y = tbl_outcomes{tbl_info(:,2)==1&tbl_info(:,3)==1,char(Y_para(y_para))};
    y = y(~(isnan(y)|isinf(y)));
    if isempty(y), continue; end

   % EF
    EFvalues = tbl_EF{tbl_info(:,2)==1&tbl_info(:,3)==1,char(calcs_x(i_calc_x))} ;
    x = EFvalues;
    mdl{1} = fitlm(x(StimSubjs),y(StimSubjs), 'linear','VarNames',[label_respondy,{label_respondx}]);
    mdl{2} =fitlm(x(ShamSubjs),y(ShamSubjs),'linear','VarNames',[label_respondy;{label_respondx}]);
 
    clear s;
    [~,~,~,~,s(1,:)] = regress(y(StimSubjs),[ones(size(x(StimSubjs))),x(StimSubjs)]);
    [~,~,~,~,s(2,:)] = regress(y(ShamSubjs),[ones(size(x(ShamSubjs))),x(ShamSubjs)]);
    try stats(1).Coefficients = mdl{1}.Coefficients; catch e, disp(['ERROR: ' e.message]); end
    try stats(2).Coefficients = mdl{2}.Coefficients; catch e, disp(['ERROR: ' e.message]); end
    try stats(1).Rsquared = mdl{1}.Rsquared; catch e, disp(['ERROR: ' e.message]); end
    try stats(2).Rsquared = mdl{2}.Rsquared; catch e, disp(['ERROR: ' e.message]); end
    stats(1).pValue = s(1,3);
    stats(2).pValue = s(2,3);
    clear s;
    close all;
    figure('Unit','centimeter','Position',[5 5 10 10],'Visible','off'); hold on; 
    x_fit = 0.9*min(x):0.005:1.1*max(x);
    ylabel( label_respondy);
    xlabel( [label_respondx '(mV/mm)']);
    for i_stim = length(stims):-1:1
        if i_stim == 1, c = colorline(i_calc_x,:); pp = StimSubjs;
        else c = [0 0 0]; pp = ShamSubjs ;end
        x_stim = x(pp);
        y_stim = y(pp);
        [p,S] = polyfit(x_stim,y_stim,1);
        y_fit = polyval(p,x_fit);
        [y_ci,y_delta] = polyconf(p,x_fit,S,'alpha',0.5);
        plot(x_fit,y_ci+y_delta,'--','LineWidth',1,'Color',c);
        plot(x_fit,y_ci-y_delta,'--','LineWidth',1,'Color',c);
        plot(x_fit,y_fit,'-','LineWidth',1,'Color',c);
        sca(i_stim) = scatter(x_stim,y_stim,12,c,...
            'Marker','o','MarkerFaceColor',c,'LineWidth',1);
    end
    set(sca(2),'MarkerFaceColor','none');
    legend(sca,'Active','Sham','Location','southwest','FontSize',6);
    xlim([min(x_fit) max(x_fit)]);

    ylim(get(gca,'ylim'));
    clear s;
    for i_stim = 1:2
        if stats(i_stim).pValue>0.05
            sign = 'n.s.';
        elseif stats(i_stim).pValue>0.01
            sign = '{\color{red}*}';
        elseif stats(i_stim).pValue>0.001
            sign = '{\color{red}**}';
        else
            sign = '{\color{red}***}';
        end     
        s{i_stim} = sprintf('%s: R^2=%.3f, adjusted R^2=%.3f, {\\itp}=%.4f (%s)', ...
            stims{i_stim},stats(i_stim).Rsquared.Ordinary,stats(i_stim).Rsquared.Adjusted,stats(i_stim).pValue,sign);
        s{i_stim} = [upper(s{i_stim}(1)), s{i_stim}(2:end)];
    end
    s = strjoin(s,'\n');
    subtitle(s,'unit','normalized','Position',[0 1.05],'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8);
    set(gca,'FontSize',8,'LineWidth',1,'Unit','centimeter','Position',[1.5 4.5 7 4]);
    fn_fig = [label_respondy label_respondx];
    fn_fig = replace(fn_fig,'\Delta','diff');
    
    if stats(1).pValue>0.05
        dn_save = [dn_results '\EFnosign'];
    else
        dn_save = [dn_results '\EFsign'];
    end
    if ~isfolder(dn_save), mkdir(dn_save); end
    print(gcf,fullfile(dn_save,[fn_fig]),'-dpdf','-r300');
    print(gcf,fullfile(dn_save,[fn_fig]),'-dpng','-r300');
    disp(['<strong>ROIs</strong>:' Sheet ])
    try disp(['#fig:' fullfile(dn_save,[fn_fig '.pdf'])]); catch e, disp('No fig'); end
    disp('===============================================================================================================');
    disp('*Real**********************************************************************************************************');
    try disp(mdl{1}); catch e, disp(['ERROR: ' e.message]); end
    fprintf('\n\n')
    disp('*Sham**********************************************************************************************************');
    try disp(mdl{2}); catch e, disp(['ERROR: ' e.message]); end            
    fprintf('\n\n\n\n\n')
        
    end
end
end
end
diary off;

%% no sperate stim and sham bold-scale
tbl_scale = readtable(fullfile(dn_scale),'VariableNamingRule','preserve');
colorline = lines(length(calcs_yy));
%% no sperate stim and sham bold-scale
% Determination of ROI
%if ~exist(fullfile(fileparts(dn_results),'stat','MLR.txt'),'file') mkdir(fullfile(fileparts(dn_results),'stat','MLR.txt'));end
%dn_results = '\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\corrBOLD\NAc_EF_abs\'
diary(fullfile(fileparts(dn_results),'stat','fMRI_Scale_individual_FC_aseg_MLR.txt'));
for y_para =  1:length(Y_para)
for yy_para = 1:length(YY_para)
for i_calc_y = 1:length(calcs_y)
for i_calc_yy = 1:length(calcs_yy)
    uyy = ''; if yy_para == 2, uyy = '%'; end
    uy = ''; if mod(y_para,2),else uy = '%'; end
    if y_para>2,  uy = ['abs',uy];end
    label_respondyy = [ '\Delta' calcs_yy{i_calc_yy} uyy ];
    label_respondyy(find(label_respondyy=='_'))='-' ;
    Sheet = char(calcs_y{i_calc_y});
    tbl_outcomes = readtable(fullfile(dn_BOLD),'Sheet',Sheet,'VariableNamingRule','preserve'); 
    label_respondy = [ '\Delta' calcs_y{i_calc_y} ')' uy ];
    label_respondy(find(label_respondy=='_'))='(&' ;
    y = tbl_outcomes{tbl_info(:,2)==1&tbl_info(:,3)==1,char(Y_para(y_para))};
    y = y(~(isnan(y)|isinf(y)));
    if isempty(y), continue; end

   % scale
    Scales = tbl_scale{tbl_info(:,2)==1&tbl_info(:,3)==1,[calcs_yy{i_calc_yy},'_2']}-tbl_scale{tbl_info(:,2)==1&tbl_info(:,3)==1,[calcs_yy{i_calc_yy},'_1']} ;
    if yy_para==2 Scales =Scales./tbl_scale{tbl_info(:,2)==1&tbl_info(:,3)==1,[calcs_yy{i_calc_yy},'_1']}; end
    yy = Scales;
    mdl{1} = fitlm(y,yy, 'linear','VarNames',[label_respondyy,{label_respondy}]);
  
 
    clear s;
    [~,~,~,~,s(1,:)] = regress(yy,[ones(size(y)),y]);
    try stats(1).Coefficients = mdl{1}.Coefficients; catch e, disp(['ERROR: ' e.message]); end  
    try stats(1).Rsquared = mdl{1}.Rsquared; catch e, disp(['ERROR: ' e.message]); end
    stats(1).pValue = s(1,3);

    clear s;
    close all;
    figure('Unit','centimeter','Position',[5 5 10 10],'Visible','off'); hold on; 
    x_fit = 0.9*min(y):0.005:1.1*max(y);
    ylabel( label_respondyy);
    xlabel( label_respondy );
    for i_stim =1% length(stims):-1:1
        if i_stim == 1, c = colorline(i_calc_yy,:); %pp = StimSubjs;
        else c = [0 0 0]; pp = ShamSubjs ;end
        x_stim = y;
        y_stim = yy;
        [p,S] = polyfit(x_stim,y_stim,1);
        y_fit = polyval(p,x_fit);
        [y_ci,y_delta] = polyconf(p,x_fit,S,'alpha',0.5);
        plot(x_fit,y_ci+y_delta,'--','LineWidth',1,'Color',c);
        plot(x_fit,y_ci-y_delta,'--','LineWidth',1,'Color',c);
        plot(x_fit,y_fit,'-','LineWidth',1,'Color',c);
        sca(i_stim) = scatter(x_stim,y_stim,12,c,...
            'Marker','o','MarkerFaceColor',c,'LineWidth',1);
    end

   % legend(sca,'Active','Sham','Location','southwest','FontSize',6);
    xlim([min(x_fit) max(x_fit)]);

    ylim(get(gca,'ylim'));
    clear s;
    for i_stim = 1%:2
        if stats(i_stim).pValue>0.05
            sign = 'n.s.';
        elseif stats(i_stim).pValue>0.01
            sign = '{\color{red}*}';
        elseif stats(i_stim).pValue>0.001
            sign = '{\color{red}**}';
        else
            sign = '{\color{red}***}';
        end     
        s{i_stim} = sprintf(':R^2=%.3f, adjusted R^2=%.3f, {\\itp}=%.4f (%s)', ...
         stats(i_stim).Rsquared.Ordinary,stats(i_stim).Rsquared.Adjusted,stats(i_stim).pValue,sign);
        s{i_stim} = [upper(s{i_stim}(1)), s{i_stim}(2:end)];
    end
    s = strjoin(s,'\n');
    subtitle(s,'unit','normalized','Position',[0 1.05],'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8);
    set(gca,'FontSize',8,'LineWidth',1,'Unit','centimeter','Position',[1.5 4.5 7 4]);
    fn_fig = [label_respondy label_respondyy];
    fn_fig = replace(fn_fig,'\Delta','diff');
    
    if stats(i_stim).pValue>0.05
        dn_save = [dn_results '\scalenosign'];
    else
        dn_save = [dn_results '\scalesign'];
    end
    if ~isfolder(dn_save), mkdir(dn_save); end
    print(gcf,fullfile(dn_save,[fn_fig]),'-dpdf','-r300');
    print(gcf,fullfile(dn_save,[fn_fig]),'-dpng','-r300');
    disp(['<strong>ROIs</strong>:' Sheet ])
    try disp(['#fig:' fullfile(dn_save,[fn_fig '.pdf'])]); catch e, disp('No fig'); end
    disp('===============================================================================================================');
        try disp(mdl{1}); catch e, disp(['ERROR: ' e.message]); end
    fprintf('\n\n')
   
        
    end
end
end
end
diary off;
