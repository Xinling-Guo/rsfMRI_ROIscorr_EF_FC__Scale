
%% load data
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
Y_para={'visit1'};
YY_para=diffratio ;
tbl_scale = readtable(fullfile(dn_scale),'VariableNamingRule','preserve');
%% no sperate stim and sham prebold-scale ： predict
RCT={'RCT','OpenL'};
cmap = [88 144 212;
    46 78 126;
    240 123 63;
    234 84 85]/255;
[ Num, Txt ] = xlsread( dn_SubjInfo,'SubjInfo' );
tbl_info =Num(:,5:8);

colorline = lines(length(calcs_yy));
ColNames = {'Stim','FC_para','FC','Sc_para','D-Scale','R2','Adjusted R2','p','Slope'};
i_stim=2;
% Determination of ROI
%if ~exist(fullfile(fileparts(dn_results),'stat','MLR.txt'),'file') mkdir(fullfile(fileparts(dn_results),'stat','MLR.txt'));end
%dn_results = '\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\corrBOLD\NAc_EF_abs\'
%Table = table([],[],[],[],[],[],'VariableNames',ColNames);
tableall = {};
k=1;
%%
for rct = 1:2
    k=1;
    Rct = RCT{rct};
    %diary(fullfile(fileparts(dn_results),'stat',['PreallFC_aseg_indi_',Rct,'.txt']));
    
    rrc=mod(rct,2);
    temp = find(tbl_info(:,2)==1&tbl_info(:,3)==1&tbl_info(:,4)==rrc);
    StimSubjs = find(tbl_info(temp,1)==1 );
    ShamSubjs = find(tbl_info(temp,1)==0 );
    tableall={};
    if rct==2 stims = {'stim'}; else stims = {'stim','sham'};end
    for y_para =  1:length(Y_para)
        for yy_para = 1:length(YY_para)
            for i_calc_y = 6%1:length(calcs_y)
                for i_calc_yy = 6%1:length(calcs_yy)
                    
                    uyy = ''; if yy_para == 2, uyy = '%'; end
                    uy = ''; if mod(y_para,2),else uy = '%'; end
                    if y_para>2,  uy = ['abs',uy];end
                    label_respondyy = [ '\Delta' calcs_yy{i_calc_yy} uyy ];
                    label_respondyy(find(label_respondyy=='_'))='-' ;
                    Sheet = char(calcs_y{i_calc_y});
                    tbl_outcomes = readtable(fullfile(dn_BOLD),'Sheet',Sheet,'VariableNamingRule','preserve');
                    label_respondy = [ 'PreAll-' calcs_y{i_calc_y} ')' uy ];
                    label_respondy(find(label_respondy=='_'))='(&' ;
                    y = tbl_outcomes{temp ,char(Y_para(y_para))};
                    y = y(~(isnan(y)|isinf(y)));
                    if isempty(y), continue; end
                    
                    % scale
                    Scales = tbl_scale{temp ,[calcs_yy{i_calc_yy},'_2']}-tbl_scale{temp ,[calcs_yy{i_calc_yy},'_1']} ;
                    if yy_para==2 Scales =Scales./tbl_scale{temp ,[calcs_yy{i_calc_yy},'_1']}; end
                    yy = Scales;
                    ind = find(isnan(yy));
                    yy(ind)=[];
                    y(ind)=[];
                    mdl= fitlm(y,yy, 'linear','VarNames',[label_respondy,{label_respondyy}]);
                    clear s;
                    [~,~,~,~,s] = regress(yy,[ones(size(y)),y]);
                    yy_z = zscore(yy);
                    y_z = zscore(y);
                    mdl_z = fitlm(y_z,yy_z, 'linear','VarNames',[label_respondy,{label_respondyy}]);
                    
                    try stats.Coefficients = mdl.Coefficients; catch e, disp(['ERROR: ' e.message]); end
                    try stats.Coefficients_z = mdl_z.Coefficients; catch e, disp(['ERROR: ' e.message]); end
                    try stats.Rsquared = mdl.Rsquared; catch e, disp(['ERROR: ' e.message]); end
                    try stats.AIC = mdl.ModelCriterion.AIC; catch, end
                    stats.pValue = s(1,3);
                                                    % cross validation
                    n_loop = 10;
                    MSE = nan(n_loop,1);
                    for i_loop = 1:n_loop, MSE(i_loop) = crossval('mse',y,yy,'Predfun',@regf,'kfold',5); end
                    stats.MSE = mean(MSE);
                    MSE =  stats(i_stim).MSE;
                    
                    clear s;
                    close all;
                    figure('Unit','centimeter','Position',[5 5 10 10],'Visible','on'); hold on;
                    x_fit = 0.9*min(y):0.005:1.1*max(y);
                    ylabel( label_respondyy);
                    xlabel( label_respondy );
                    
                    c = colorline(i_calc_yy,:);
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
                    %legend(sca,'Active','Sham','Location','southwest','FontSize',6);
                    xlim([min(x_stim)-0.02*abs(min(x_stim)) max(x_stim)+0.02*abs(max(x_stim))]);
                    
                    
                    ylim(get(gca,'ylim'));
                    clear s;
                    beta = stats.Coefficients_z{2,'Estimate'};
                    pValue = stats.Coefficients_z{2,'pValue'};
                    s = sprintf('R^2=%.2f, Adj. R^2=%.2f, MSE=%.1f\nAIC=%.2f, {\\itp}=%.3f (%s), {\\beta}=%.2f, {\\itp}=%.3f (%s)', stats.Rsquared.Ordinary, stats.Rsquared.Adjusted,stats.MSE,stats.AIC,stats.pValue, sign(stats.pValue),beta,pValue,sign(pValue));
                    
                    
                    subtitle(s,'unit','normalized','Position',[0 1.05],'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8);
                    set(gca,'FontSize',8,'LineWidth',1,'Unit','centimeter','Position',[1.5 4.5 7 4]);
                    fn_fig = [Rct label_respondy label_respondyy];
                    fn_fig = replace(fn_fig,'\Delta','diff');
                    
                    
                    if stats.pValue>0.05
                        dn_save = [dn_results '\Preallnosign'];
                    else
                        dn_save = [dn_results '\Preallsign'];
                        i_stim=1;
                        tableall(k,:) ={i_stim,Y_para{y_para},label_respondy,YY_para{yy_para}, label_respondyy,stats.Rsquared.Ordinary,...
                            stats.Rsquared.Adjusted,stats.pValue,table2array(mdl.Coefficients(2,1))};
                        
                        k=k+1;
                    end
                    if ~isfolder(dn_save), mkdir(dn_save); end
                    print(gcf,fullfile(dn_save,[fn_fig]),'-dpdf','-r300');
                    %print(gcf,fullfile(dn_save,[fn_fig]),'-dpng','-r300');
                    disp(['<strong>ROIs</strong>:' calcs_y{i_calc_y} ])
                    try disp(['#fig:' fullfile(dn_save,[fn_fig '.pdf'])]); catch e, disp('No fig'); end
                    disp('===============================================================================================================');
                    disp('*All**********************************************************************************************************');
                    try disp(mdl); catch e, disp(['ERROR: ' e.message]); end
                    fprintf('\n\n')
                    fprintf('\n\n\n\n\n')
                    
                end
            end
        end
    end
    diary off;
    if length(tableall)>0
        Table=cell2table(tableall,'VariableNames',ColNames);
        tablefile = [dn_results,'\onePDF\','presign.xls'];
        writetable(Table,tablefile,'Sheet',[Rct,'PreallFC_scale']);
    end
end
%% sperate stim and sham prebold-scale ： predict
colorline = lines(length(calcs_yy));
ColNames = {'Stim','FC_para','FC','Sc_para','D-Scale','R2','Adjusted R2','p','Slope'};
for rct = 2%1:2
    Rct = RCT{rct};
    diary(fullfile(fileparts(dn_results),'stat',['PreFC_aseg_indi_',Rct,'.txt']));
    rrc=mod(rct,2);
    temp = find(tbl_info(:,2)==1&tbl_info(:,3)==1&tbl_info(:,4)==rrc);
    StimSubjs = find(tbl_info(temp,1)==1 );
    ShamSubjs = find(tbl_info(temp,1)==0 );
    if rct==2 stims = {'stim'}; else stims = {'stim','sham'};end
    tableall={};
    k=1;
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
                    label_respondy = [ 'Pre-' calcs_y{i_calc_y} ')' uy ];
                    label_respondy(find(label_respondy=='_'))='(&' ;
                    y = tbl_outcomes{temp,char(Y_para(y_para))};
                    y = y(~(isnan(y)|isinf(y)));
                    if isempty(y), continue; end
                    
                    % scale
                    Scales = tbl_scale{temp,[calcs_yy{i_calc_yy},'_2']}-tbl_scale{temp,[calcs_yy{i_calc_yy},'_1']} ;
                    if yy_para==2 Scales =Scales./tbl_scale{temp,[calcs_yy{i_calc_yy},'_1']}; end
                    yy = Scales;
                    
                    figure('Unit','centimeter','Position',[5 5 20 10],'Visible','off'); hold on;
                    x_fit = 0.9*min(y):0.005:1.1*max(y);
                    for i_stim = length(stims):-1:1
                        if i_stim == 1, c = colorline(i_calc_yy,:); pp = StimSubjs;
                        else c = [0 0 0]; pp = ShamSubjs ;end
                        subplot(1,2,(i_stim));
                        x_stim = y(pp);
                        y_stim = yy(pp);
                        ind = find(isnan(y_stim));
                        x_stim(ind)=[];
                        y_stim(ind)=[];
                        %ismember(ind,StimSubjs)
                        mdl{i_stim} = fitlm(x_stim,y_stim, 'linear','VarNames',[label_respondy,{label_respondyy}]);
                        clear s;
                        [~,~,~,~,s(i_stim,:)] = regress(y_stim,[ones(size(x_stim)),x_stim]);
                        x_stim_z = zscore( x_stim);
                        y_stim_z = zscore( y_stim);
                        mdl_z{i_stim} = fitlm(x_stim_z,y_stim_z , 'linear','VarNames',[label_respondy,{label_respondyy}]);
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
                        [p,S] = polyfit(x_stim,y_stim,1);
                        y_fit = polyval(p,x_fit);
                        [y_ci,y_delta] = polyconf(p,x_fit,S,'alpha',0.5);
                        plot(x_fit,y_ci+y_delta,'--','LineWidth',1,'Color',c); hold on;
                        plot(x_fit,y_ci-y_delta,'--','LineWidth',1,'Color',c);; hold on;
                        plot(x_fit,y_fit,'-','LineWidth',1,'Color',c); hold on;
                        sca(i_stim) = scatter(x_stim,y_stim,12,c,...
                            'Marker','o','MarkerFaceColor',c,'LineWidth',1); hold on;
                        xlim([min(x_stim)-0.02*abs(min(x_stim)) max(x_stim)+0.02*abs(max(x_stim))]);
                        ylim(get(gca,'ylim'));
                        ylabel( label_respondyy);
                        xlabel( label_respondy );
                        clear s;
                        beta = stats(i_stim).Coefficients_z{2,'Estimate'};
                        pValue = stats(i_stim).Coefficients_z{2,'pValue'};
                        s = sprintf('R^2=%.2f, Adj. R^2=%.2f, MSE=%.1f\nAIC=%.2f, {\\itp}=%.3f (%s), {\\beta}=%.2f, {\\itp}=%.3f (%s)', stats(i_stim).Rsquared.Ordinary, stats(i_stim).Rsquared.Adjusted,stats(i_stim).MSE,stats(i_stim).AIC,stats(i_stim).pValue, sign(stats(i_stim).pValue),beta,pValue,sign(pValue));
                        subtitle(s,'unit','normalized','Position',[0 1.001],'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8);
                    end
                    fn_fig = [Rct label_respondy label_respondyy];
                    fn_fig = replace(fn_fig,'\Delta','diff');
                    
                    if stats(1).pValue>0.05
                        dn_save = [dn_results '\PreFCnosign'];
                    else
                        dn_save = [dn_results '\PreFCsign'];
                        i_stim=1;
                        tableall(k,:) ={i_stim,Y_para{y_para},label_respondy,YY_para{yy_para}, label_respondyy,stats(1).Rsquared.Ordinary,...
                            stats(1).Rsquared.Adjusted,stats(1).pValue,table2array(mdl{1}.Coefficients(2,1))};
                        
                        k=k+1;
                    end
                    if ~isfolder(dn_save), mkdir(dn_save); end
                    print(gcf,fullfile(dn_save,[fn_fig]),'-dpdf','-r300');
                    %print(gcf,fullfile(dn_save,[fn_fig]),'-dpng','-r300');
                    disp(['<strong>ROIs</strong>:' calcs_y{i_calc_y} ])
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
    Table=cell2table(tableall,'VariableNames',ColNames);
    tablefile =  [dn_results,'\onePDF\','presign.xls'];
    writetable(Table,tablefile,'Sheet',[Rct,'_PreFC_scale']);
end
%% sperate stim and sham prebold+EF - scale ： predict（save sign to excle）
% Determination of ROI
%if ~exist(fullfile(fileparts(dn_results),'stat','MLR.txt'),'file') mkdir(fullfile(fileparts(dn_results),'stat','MLR.txt'));end
RCT={'RCT','OpenL'};
cmap = [88 144 212;
    46 78 126;
    240 123 63;
    234 84 85]/255;
[ Num, Txt ] = xlsread( dn_SubjInfo,'SubjInfo' );
tbl_info =Num(:,5:8);
colorline = lines(length(calcs_yy));
ColNames = {'Stim','EF_para','EF','FC_para','preFC','Sc_para','D-Scale','R2','Adjusted R2','p','x1-Slope','x2-Slope',...
   'MSE', 'x1-beta','x1-p','x2-beta','x2-p'};
for rct = 1:2
    Rct = RCT{rct};
    diary(fullfile(fileparts(dn_results),'stat',['PreEFFC_aseg_indi_',Rct,'.txt']));
    rrc=mod(rct,2);
    temp = find(tbl_info(:,2)==1&tbl_info(:,3)==1&tbl_info(:,4)==rrc);
    StimSubjs = find(tbl_info(temp,1)==1 );
    ShamSubjs = find(tbl_info(temp,1)==0 );
    if rct==2 stims = {'stim'}; else stims = {'stim','sham'};end
    tableall={};
    k=1;
    for x_para =  1:length(X_para)
        tbl_EF = readtable(fullfile(dn_EF),'Sheet',X_para{x_para},'VariableNamingRule','preserve');
        for y_para =  1:length(Y_para)
            for yy_para = 1:length(YY_para)
                for i_calc_x = 1:length(calcs_x)
                    for i_calc_y = 1:length(calcs_y)
                        for i_calc_yy = 1:length(calcs_yy)
                            uyy = ''; if yy_para == 2, uyy = '%'; end
                            uy = ''; if mod(y_para,2),else uy = '%'; end
                            if y_para>2,  uy = ['abs',uy];end
                            label_respondx = [ calcs_x{i_calc_x} '-' X_para{x_para} ];
                            label_respondx(find(label_respondx=='_'))='-' ;
                            label_respondyy = [ '\Delta' calcs_yy{i_calc_yy} uyy ];
                            label_respondyy(find(label_respondyy=='_'))='-' ;
                            Sheet = char(calcs_y{i_calc_y});
                            tbl_outcomes = readtable(fullfile(dn_BOLD),'Sheet',Sheet,'VariableNamingRule','preserve');
                            label_respondy = [ 'Pre-' calcs_y{i_calc_y} ')' uy ];
                            label_respondy(find(label_respondy=='_'))='(&' ;
                            y = tbl_outcomes{temp,char(Y_para(y_para))};
                            y = y(~(isnan(y)|isinf(y)));
                            FC = y;
                            if isempty(y), continue; end
                            % EF
                            EFvalues = tbl_EF{temp,char(calcs_x(i_calc_x))} ;
                            x = EFvalues;
                            EF =x ;
                            % scale
                            Scales = tbl_scale{temp,[calcs_yy{i_calc_yy},'_2']}-tbl_scale{temp,[calcs_yy{i_calc_yy},'_1']} ;
                            if yy_para==2 Scales =Scales./tbl_scale{temp,[calcs_yy{i_calc_yy},'_1']}; end
                            yy = Scales;
                            % X
                            load carsmall
                            X = [EF,FC];
                            figure('Unit','centimeter','Position',[5 5 20 10],'Visible','off'); hold on;
                            % plot
                            for i_stim = length(stims):-1:1
                                if i_stim == 1, c = colorline(i_calc_yy,:); pp = StimSubjs;
                                else c = [0 0 0]; pp = ShamSubjs ;end
                                subplot(1,2,(i_stim));
                                x_stim = EF(pp);
                                y_stim = FC(pp);
                                z_stim = Scales(pp);
                                ind = find(isnan(x_stim));
                                ind = [ind find(isnan(z_stim))];
                                x_stim(ind)=[];
                                y_stim(ind)=[];
                                z_stim(ind)=[];
                                X_stim = [x_stim y_stim];
                                x1=x_stim;
                                x2=y_stim;
                                mdl{i_stim} = fitlm(X_stim,z_stim, 'linear','VarNames',[label_respondx,label_respondy,{label_respondyy}]);
                                clear s;
                                X2_stim = [ones(size(x1)) x1 x2 x1.*x2];
                                [b,BINT,R,RINT,s(i_stim,:)] = regress(z_stim,X2_stim);
                                X_stim_z = zscore( X_stim);
                                z_stim_z = zscore( z_stim);
                                mdl_z{i_stim} = fitlm(X_stim_z,z_stim_z , 'linear','VarNames',[label_respondx,label_respondy,{label_respondyy}]);
                                try stats(i_stim).Coefficients = mdl{i_stim}.Coefficients; catch e, disp(['ERROR: ' e.message]); end
                                try stats(i_stim).Coefficients_z = mdl_z{i_stim}.Coefficients; catch e, disp(['ERROR: ' e.message]); end
                                try stats(i_stim).Rsquared = mdl{i_stim}.Rsquared; catch e, disp(['ERROR: ' e.message]); end
                                try stats(i_stim).AIC = mdl{i_stim}.ModelCriterion.AIC; catch, end
                                stats(i_stim).pValue = s(i_stim,3);
                                % cross validation
                                % cross validation
                                n_loop = 10;
                                MSE = nan(n_loop,1);
                                for i_loop = 1:n_loop, MSE(i_loop) = crossval('mse',X_stim,z_stim,'Predfun',@regf,'kfold',5); end
                                stats(i_stim).MSE = mean(MSE);
                                MSE =  stats(i_stim).MSE;
                                clear s;
                                scatter3(x1,x2,z_stim ,'filled','MarkerFaceColor',c);
                                hold on  %在刚刚那副散点图上接着画
                                % if i_stim == 1
                                x1fit = 0.9*min(x1):0.005:1.1*max(x1);;   %设置x1的数据间隔
                                x2fit = 0.9*min(x2):0.005:1.1*max(x2);     %设置x2的数据间隔
                                [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);    %生成一个二维网格平面，也可以说生成X1FIT,X2FIT的坐标
                                hold on;
                                YFIT=b(1)+b(2)*X1FIT+b(3)*X2FIT+b(4)*X1FIT.*X2FIT;    %代入已经求得的参数，拟合函数式
                                mesh(X1FIT,X2FIT,YFIT);
                                ylabel( label_respondy);
                                xlabel( label_respondx );
                                zlabel(label_respondyy);
                                %set(get(gca,'xlabel'),'rotation',25,'Position',[0.04 -0.65 -5]);
                                %set(get(gca,'ylabel'),'rotation',-35,'Position',[0.0004 0.5 -5]);
                                h = rotate3d;
                                set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
                                set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
                                set(gcf, 'ResizeFcn', @align_axislabel)
                                align_axislabel([], gca)
                                axislabel_translation_slider;
                                % end
                                clear s;
                                beta1 = stats(i_stim).Coefficients_z{2,'Estimate'};
                                pValue1 = stats(i_stim).Coefficients_z{2,'pValue'};
                                beta2 = stats(i_stim).Coefficients_z{3,'Estimate'};
                                pValue2 = stats(i_stim).Coefficients_z{3,'pValue'};
                                s = sprintf('R^2=%.2f, Adj. R^2=%.2f, MSE=%.1f, AIC=%.2f, {\\itp}=%.3f (%s)\nx1:\\beta=%.2f, {\\itp}=%.3f (%s),x2:\\beta=%.2f, {\\itp}=%.3f (%s)',...
                                    stats(i_stim).Rsquared.Ordinary, stats(i_stim).Rsquared.Adjusted,stats(i_stim).MSE,stats(i_stim).AIC,stats(i_stim).pValue, sign(stats(i_stim).pValue),...
                                    beta1,pValue1,sign(pValue1),beta2,pValue2,sign(pValue2) );
                                subtitle(s,'unit','normalized','Position',[0 1.001],'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8);
                            end
                            fn_fig = [Rct label_respondx label_respondy label_respondyy];
                            fn_fig = replace(fn_fig,'\Delta','diff');
                            
                            if stats(1).pValue>0.05
                                dn_save = [dn_results '\PreEFFCnosign'];
                            else
                                dn_save = [dn_results '\PreEFFCsign'];
                                i_stim=1;
                                tableall(k,:) ={i_stim,X_para{x_para},label_respondx,Y_para{y_para},label_respondy,YY_para{yy_para}, label_respondyy,stats(1).Rsquared.Ordinary,...
                                    stats(1).Rsquared.Adjusted,stats(1).pValue,table2array(mdl{1}.Coefficients(2,1)),table2array(mdl{1}.Coefficients(3,1)),...
                                   stats(i_stim).MSE,beta1,pValue1,beta2,pValue2 };
                                k=k+1;
                                
                            end
                            if ~isfolder(dn_save), mkdir(dn_save); end
                            print(gcf,fullfile(dn_save,[fn_fig]),'-dpdf','-r300');
                            %print(gcf,fullfile(dn_save,[fn_fig]),'-dpng','-r300');
                            disp(['<strong>ROIs</strong>:' calcs_y{i_calc_y} ])
                            try disp(['#fig:' fullfile(dn_save,[fn_fig '.pdf'])]); catch e, disp('No fig'); end
                            disp('===============================================================================================================');
                            disp('*Real**********************************************************************************************************');
                            try disp(mdl{1}); catch e, disp(['ERROR: ' e.message]); end
                            fprintf('\n\n')
                            disp('*Sham**********************************************************************************************************');
                            try disp(mdl{2}); catch e, disp(['ERROR: ' e.message]); end
                            fprintf('\n\n\n\n\n')
                            all_fig = findall(0, 'type', 'figure');
                            close(all_fig)
                        end
                    end
                end
            end
        end
    end
    
    diary off;
    Table=cell2table(tableall,'VariableNames',ColNames);
    tablefile =  [dn_results,'\onePDF\','presign.xls'];
    writetable(Table,tablefile,'Sheet',[Rct,'_PreEFFC_scale']);
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