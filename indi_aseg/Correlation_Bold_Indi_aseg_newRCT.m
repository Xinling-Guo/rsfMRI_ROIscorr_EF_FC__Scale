%%
imgDir=dir('*.*');
sheet=2;
[ Num, Txt ] = xlsread( 'D:\jingzongfmri\stimTFnew.xlsx',sheet );
filename = Txt(2:46,[3,4]);
subjnum = length(filename);
Visit = {'v1','v2'};
%% FC excle
% excle for pre and post FC
%FCfile='\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\tbl\fMRI_individual_FC_aseg.xlsx';
FCfile='\\172.16.4.19\eeg\Projects\tDCS\202201-SMHC_cathode_FP2\tbl\fMRI_individual_FC_aseg_RCT_new.xlsx';

%% set session 

ROIs = {'2012','2014','1027','2027','1028','2028','26','58','1008','2008'};
ROIsName ={'rl-OFC','rm-OFC','ll-PFC','rl-PFC','lm-PFC','rm-PFC','l-NAc','r-NAc','l-IPL','r-IPL'};
FCpath1 ='D:\jingzongfmri\data\fmri_indi_';
FCpath2 ='_new\Results\ROISignals_FunImgARCFS\';
%la:lateral 外侧；me:medial 内侧 IPL :inferiorparietal
%%  write FC （post-pre） to excle
for roi1 = 1:length(ROIs)
    for roi2 = roi1+1:length(ROIs)
        Sheetname = ['FC_',char(ROIsName(roi1)),'_',char(ROIsName(roi2))];
        ColNames = {'index','visit1','visit2','D-value','D-ratio','absD-value','absD-ratio'};
        Table  = table();
        FCRoisAllSubj = zeros((subjnum),length(Visit));
        for subj = 40:subjnum
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
        absdifferent = absornot(FCRoisAllSubj(:,2))-absornot(FCRoisAllSubj(:,1));
        differentR = different./FCRoisAllSubj(:,1);
        absdifferentR = absdifferent./absornot(FCRoisAllSubj(:,1));
        Table=table([1:subjnum]',FCRoisAllSubj(:,1),FCRoisAllSubj(:,2),different,differentR,absdifferent,absdifferentR,'VariableNames',ColNames);
        writetable(Table(40:45,:),FCfile,'Sheet',Sheetname,'WriteMode','Append');
    end
end



%% Multivarables linear regression for BOLD-ef

clear;
clc;
close all;

dn_results = '\\172.16.4.19\eeg\Projects\tDCS\202201-SMHC_cathode_FP2\data\fMRI\CorrBold_individual_aseg_new\'
dn_root = '\\172.16.4.19\eeg\Projects\tDCS\202201-SMHC_cathode_FP2\tbl\';
dn_EF = [dn_root,'EF.xlsx'];
dn_BOLD =  [dn_root,'fMRI_individual_FC_aseg_RCT_new.xlsx'];
dn_scale = [dn_root,'scales.xlsx'];
dn_SubjInfo =  [dn_root,'fMRI_info.xlsx'];
[Type,Sheet,Format]=xlsfinfo(dn_BOLD);
calcs_x = {'rl-OFC','rm-OFC','rl-PFC','rm-PFC'};
calcs_y= Sheet;
calcs_yy = {'PSQI','YBOCS','YBOCS_Obsessions','YBOCS_Compulsions','BDI','BAI'};

diffratio = {'D-value','D-ratio'};
absornot={'','abs'};
stims = {'stim','sham'};
labels_stim = [1 0];
colorline = lines(length(calcs_x));
X_para={'mean_EF','perc95_EF'};
Y_para={'D-value',	'D-ratio',	'absD-value',	'absD-ratio'};
YY_para=diffratio ;
%% sperate stim and sham EF-bold RCT or OpenL
RCT={'RCT','OpenL'};
cmap = [88 144 212;
        46 78 126;
        240 123 63;
        234 84 85]/255;
[ Num, Txt ] = xlsread( dn_SubjInfo,'SubjInfo' );
tbl_info =Num(:,5:8);
ColNames = {'Stim','EF_para','EF','FC_para','D-FC','R2','Adjusted R2','p','Slope','MSE'};

% Determination of ROI
%if ~exist(fullfile(fileparts(dn_results),'stat','MLR.txt'),'file') mkdir(fullfile(fileparts(dn_results),'stat','MLR.txt'));end
for rct = 2%1:2
    k=1;
    Rct = RCT{rct};
    diary(fullfile(fileparts(dn_results),'stat',['EF_FC_aseg_indi_',Rct,'.txt']));

    rrc=mod(rct,2);
    temp = find(tbl_info(:,2)==1&tbl_info(:,3)==1&tbl_info(:,4)==rrc);
    StimSubjs = find(tbl_info(temp,1)==1 );
    ShamSubjs = find(tbl_info(temp,1)==0 );
    if rct==2 stims = {'stim'}; else stims = {'stim','sham'};end
    tableall={};
    for x_para =  1:length(X_para)
        tbl_EF = readtable(fullfile(dn_EF),'Sheet',X_para{x_para},'VariableNamingRule','preserve');
        tbl_EF([39,40],:)=[];
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
        FC = tbl_outcomes{temp,char(Y_para(y_para))};
        FC = FC(~(isnan(FC)|isinf(FC)));
        if isempty(FC), continue; end

       % EF
        EF = tbl_EF{temp,char(calcs_x(i_calc_x))} ;
        figure('Unit','centimeter','Position',[5 5 21 12],'Visible','off'); hold on; 
        for i_stim = length(stims):-1:1
            if i_stim == 1, c = colorline(i_calc_x,:); pp = StimSubjs;
            else c = [0 0 0]; pp = ShamSubjs ;end
            x_stim = EF(pp);
            y_stim = FC(pp);
            ind = find(isnan(x_stim));
            x_stim(ind)=[];
                               y_stim(ind)=[];
            
            mdl{i_stim} = fitlm(x_stim,y_stim, 'linear','VarNames',[label_respondx,{label_respondy}]);
            clear s;
            [~,~,~,~,s(i_stim,:)] = regress(y_stim,[ones(size(x_stim)),x_stim]);
            % calculate the β value of ROIs
            y_z = zscore(y_stim);
            x_z = zscore(x_stim);
            mdl_z{i_stim} = fitlm(x_z,y_z, 'linear','VarNames',[label_respondx,{label_respondy}]);
             % summary the statistic values   
            try stats(i_stim).Coefficients = mdl{i_stim}.Coefficients; catch e, disp(['ERROR: ' e.message]); end
            try stats(i_stim).Coefficients_z = mdl_z{i_stim}.Coefficients; catch e, disp(['ERROR: ' e.message]); end
            try stats(i_stim).Rsquared = mdl{i_stim}.Rsquared; catch e, disp(['ERROR: ' e.message]); end
            try stats(i_stim).AIC = mdl{i_stim}.ModelCriterion.AIC; catch, end
            stats(i_stim).pValue = s(i_stim,3);
            % cross validation
            n_loop = 10;
            MSE = nan(n_loop,1);
            for i_loop = 1:n_loop, MSE(i_loop) = crossval('mse',x_stim,y_stim,'Predfun',@regf,'kfold',5); end
            stats(i_stim).MSE = mean(MSE);
            MSE =  stats(i_stim).MSE;
            clear s;
            subplot(1,2,(i_stim));
            [p,S] = polyfit(x_stim,y_stim,1);
            x_fit = 0.9*min(x_stim):0.005:1.5*max(x_stim);
            y_fit = polyval(p,x_fit);
            [y_ci,y_delta] = polyconf(p,x_fit,S,'alpha',0.5);
            plot(x_fit,y_ci+y_delta,'--','LineWidth',1,'Color',c);hold on;
            plot(x_fit,y_ci-y_delta,'--','LineWidth',1,'Color',c); hold on;
            plot(x_fit,y_fit,'-','LineWidth',1,'Color',c); hold on;
            sca(i_stim) = scatter(x_stim,y_stim,12,c,...
                'Marker','o','MarkerFaceColor',c,'LineWidth',1);
            xlim([min(x_stim)-0.02*abs(min(x_stim)) max(x_stim)+0.02*abs(max(x_stim))]);

            ylabel( label_respondy);
            xlabel( [label_respondx '(mV/mm)']);
            legend off;
            title('');
            beta = stats(i_stim).Coefficients_z{2,'Estimate'};
            pValue = stats(i_stim).Coefficients_z{2,'pValue'};
            s = sprintf('R^2=%.2f, Adj. R^2=%.2f, MSE=%.1f\nAIC=%.2f, {\\itp}=%.3f (%s), {\\beta}=%.2f, {\\itp}=%.3f (%s)', stats(i_stim).Rsquared.Ordinary, stats(i_stim).Rsquared.Adjusted,stats(i_stim).MSE,stats(i_stim).AIC,stats(i_stim).pValue, sign(stats(i_stim).pValue),beta,pValue,sign(pValue));
            subtitle(s,'Unit','normalize','Position',[0 1.01], ...
                'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8); 
        end
        fn_fig = [Rct label_respondy label_respondx];
        fn_fig = replace(fn_fig,'\Delta','diff');
        if stats(1).pValue>0.05
            dn_save = [dn_results '\EF','nosign'];
        else
        dn_save = [dn_results '\EF','sign'];
        i_stim=1;
        tableall(k,:) ={i_stim,X_para{x_para},label_respondx,Y_para{y_para},label_respondy,stats(1).Rsquared.Ordinary,...
            stats(1).Rsquared.Adjusted,stats(1).Coefficients_z{2,'pValue'},table2array(mdl{1}.Coefficients(2,1)),MSE};
        
        k=k+1;
        end
        if ~isfolder(dn_save), mkdir(dn_save); end
        print(gcf,fullfile(dn_save,[fn_fig]),'-dpdf','-r300');
        %print(gcf,fullfile(dn_save,[fn_fig]),'-dpng','-r300');
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
    if length(tableall)>0
        Table=cell2table(tableall,'VariableNames',ColNames);
        tablefile = [dn_results,'\onePDF\','EF_FC_Scale_sign.xls'];
        writetable(Table,tablefile,'Sheet',[Rct,'EF_FC_sign']);
    end
end
%% no sperate stim and sham bold-scale
tbl_scale = readtable(fullfile(dn_scale),'VariableNamingRule','preserve');
colorline = lines(length(calcs_yy));
tbl_scale([39,40],:)=[];
%% no sperate stim and sham bold-scale
% Determination of ROI
%if ~exist(fullfile(fileparts(dn_results),'stat','MLR.txt'),'file') mkdir(fullfile(fileparts(dn_results),'stat','MLR.txt'));end
%dn_results = '\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\corrBOLD\NAc_EF_abs\'
RCT={'RCT','OpenL'};
cmap = [88 144 212;
        46 78 126;
        240 123 63;
        234 84 85
        70  120 200
        140  220 100]/255;
[ Num, Txt ] = xlsread( dn_SubjInfo,'SubjInfo' );
tbl_info =Num(:,5:8);
ColNames = {'Stim','FC_para','FC','Sc_para','D-Sc','R2','Adjusted R2','p','Slope','MSE'};

for rct = 2%1:2
    k=1;
    Rct = RCT{rct};
    diary(fullfile(fileparts(dn_results),'stat',[Rct,'FC_Scale_aseg_indi_RCT.txt']));
    rrc=mod(rct,2);
    temp = find(tbl_info(:,2)==1&tbl_info(:,3)==1&tbl_info(:,4)==rrc);
    tableall={};
    if rct==2 stims = {'stim'}; else stims = {'stim','sham'};end
    for y_para =  1:length(Y_para)
    for yy_para = 1:length(YY_para)
    for i_calc_y = 1:length(calcs_y)
    for i_calc_yy = 1:length(calcs_yy)
        istim = 1;
        uyy = ''; if yy_para == 2, uyy = '%'; end
        uy = ''; if mod(y_para,2),else uy = '%'; end
        if y_para>2,  uy = ['abs',uy];end
        label_respondyy = [ '\Delta' calcs_yy{i_calc_yy} uyy ];
        label_respondyy(find(label_respondyy=='_'))='-' ;
        Sheet = char(calcs_y{i_calc_y});
        tbl_outcomes = readtable(fullfile(dn_BOLD),'Sheet',Sheet,'VariableNamingRule','preserve'); 
        label_respondy = [ '\Delta' calcs_y{i_calc_y} ')' uy ];
        label_respondy(find(label_respondy=='_'))='(&' ;
        y_stim = tbl_outcomes{temp,char(Y_para(y_para))};
        y_stim = y_stim(~(isnan(y_stim)|isinf(y_stim)));
        if isempty(y_stim), continue; end

       % scale
        Scales = tbl_scale{temp,[calcs_yy{i_calc_yy},'_2']}-tbl_scale{temp,[calcs_yy{i_calc_yy},'_1']} ;
        if yy_para==2 Scales =Scales./tbl_scale{temp,[calcs_yy{i_calc_yy},'_1']}; end
        yy = Scales;
        ind = find(isnan(yy));
        yy(ind)=[];
        y_stim(ind)=[];
        mdl{1} = fitlm(y_stim,yy, 'linear','VarNames',[label_respondy,{label_respondyy}]);


        clear s;
        [~,~,~,~,s(1,:)] = regress(yy,[ones(size(y_stim)),y_stim]);
         % calculate the β value of ROIs
        yy_z = zscore(yy);
        y_z = zscore(y_stim);
        mdl_z{1} = fitlm(y_z,yy_z, 'linear','VarNames',[label_respondy,{label_respondyy}]);

        try stats(1).Coefficients = mdl{1}.Coefficients; catch e, disp(['ERROR: ' e.message]); end  
        try stats(1).Coefficients_z = mdl_z{1}.Coefficients; catch e, disp(['ERROR: ' e.message]); end  
        try stats(1).Rsquared = mdl{1}.Rsquared; catch e, disp(['ERROR: ' e.message]); end
        try stats(1).AIC = mdl{1}.ModelCriterion.AIC; catch, end   
        stats(1).pValue = s(1,3);
        % cross validation
                    n_loop = 10;
            MSE = nan(n_loop,1);
            for i_loop = 1:n_loop, MSE(i_loop) = crossval('mse',y_stim,yy,'Predfun',@regf,'kfold',5); end
            stats(1).MSE = mean(MSE);
            MSE =  stats(1).MSE;

        clear s;
        close all;
        figure('Unit','centimeter','Position',[5 5 10 10],'Visible','off'); hold on; 
        x_fit = 0.9*min(y_stim):0.005:1.1*max(y_stim);
        ylabel( label_respondyy);
        xlabel( label_respondy );
        i_stim =1;
            if i_stim == 1, c = colorline(i_calc_yy,:); %pp = StimSubjs;
            else c = [0 0 0]; pp = ShamSubjs ;end
            x_stim = y_stim;
            y_stim = yy;
            [p,S] = polyfit(x_stim,y_stim,1);
            y_fit = polyval(p,x_fit);
            [y_ci,y_delta] = polyconf(p,x_fit,S,'alpha',0.5);
            plot(x_fit,y_ci+y_delta,'--','LineWidth',1,'Color',c);
            plot(x_fit,y_ci-y_delta,'--','LineWidth',1,'Color',c);
            plot(x_fit,y_fit,'-','LineWidth',1,'Color',c);
            sca(i_stim) = scatter(x_stim,y_stim,12,c,...
                'Marker','o','MarkerFaceColor',c,'LineWidth',1);
        

       % legend(sca,'Active','Sham','Location','southwest','FontSize',6);
        xlim([min(x_stim)-0.02*abs(min(x_stim)) max(x_stim)+0.02*abs(max(x_stim))]);

        ylim(get(gca,'ylim'));
        clear s;
        beta = stats(i_stim).Coefficients_z{2,'Estimate'};
            pValue = stats(i_stim).Coefficients_z{2,'pValue'};
            s = sprintf('R^2=%.2f, Adj. R^2=%.2f, MSE=%.1f\nAIC=%.2f, {\\itp}=%.3f (%s), {\\beta}=%.2f, {\\itp}=%.3f (%s),MSE=%.3f', stats(i_stim).Rsquared.Ordinary, stats(i_stim).Rsquared.Adjusted,stats(i_stim).MSE,stats(i_stim).AIC,stats(i_stim).pValue, sign(stats(i_stim).pValue),beta,pValue,sign(pValue),MSE);
     
        %s = strjoin(s,'\n');
        subtitle(s,'unit','normalized','Position',[0 1.05],'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8);
        set(gca,'FontSize',8,'LineWidth',1,'Unit','centimeter','Position',[1.5 4.5 7 4]);
        fn_fig = [Rct label_respondy label_respondyy];
        fn_fig = replace(fn_fig,'\Delta','diff');

        if stats(i_stim).pValue>0.05
            dn_save = [dn_results '\scalenosign'];
        else
            dn_save = [dn_results '\scalesign'];
           i_stim=1;
        tableall(k,:) ={i_stim,Y_para{y_para},label_respondy,YY_para{yy_para},label_respondyy,stats(1).Rsquared.Ordinary,...
            stats(1).Rsquared.Adjusted,stats(1).Coefficients_z{2,'pValue'},table2array(mdl{1}.Coefficients(2,1)),MSE};
        
        k=k+1;
        end
        if ~isfolder(dn_save), mkdir(dn_save); end
        print(gcf,fullfile(dn_save,[fn_fig]),'-dpdf','-r300');
%         print(gcf,fullfile(dn_save,[fn_fig]),'-dpng','-r300');
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
    if length(tableall)>0
        Table=cell2table(tableall,'VariableNames',ColNames);
        tablefile = [dn_results,'\onePDF\','EF_FC_Scale_sign.xls'];
        writetable(Table,tablefile,'Sheet',[Rct,'_FC_Scale_sign']);
    end
end
%%
function yfit = regf(Xtrain,ytrain,Xtest)
mdl = fitlm(Xtrain,ytrain);
yfit = predict(mdl,Xtest);
end

function s = sign(p)
    if p>0.05
        s = 'n.s.';
    elseif p>0.01
        s = '{\color{red}*}';
    elseif p>0.001
        s = '{\color{red}**}';
    else
        s = '{\color{red}***}';
    end     
end