
%% load data
clear;
clc;
close all;

dn_results = '\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\CorrBold_individual_aseg\';
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
absornot={'','abs'};
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
Y2_para=absornot;
Y_para={'visit1'};
YY_para=diffratio ;
tbl_scale = readtable(fullfile(dn_scale),'VariableNamingRule','preserve');
%% no sperate stim and sham prebold-scale ： predict
colorline = lines(length(calcs_yy));
ColNames = {'Stim','FC_para','FC','Sc_para','D-Scale','R2','Adjusted R2','p','Slope'};
i_stim=2;
% Determination of ROI
%if ~exist(fullfile(fileparts(dn_results),'stat','MLR.txt'),'file') mkdir(fullfile(fileparts(dn_results),'stat','MLR.txt'));end
%dn_results = '\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\corrBOLD\NAc_EF_abs\'
diary(fullfile(fileparts(dn_results),'stat','PreAll_individual_FC_aseg.txt'));
%Table = table([],[],[],[],[],[],'VariableNames',ColNames);
tableall = {};
k=1;
%%
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
    label_respondy = [ 'PreAll-' calcs_y{i_calc_y} ')' uy ];
    label_respondy(find(label_respondy=='_'))='(&' ;
    y = tbl_outcomes{tbl_info(:,2)==1&tbl_info(:,3)==1,char(Y_para(y_para))};
    y = y(~(isnan(y)|isinf(y)));
    if isempty(y), continue; end

   % scale
    Scales = tbl_scale{tbl_info(:,2)==1&tbl_info(:,3)==1,[calcs_yy{i_calc_yy},'_2']}-tbl_scale{tbl_info(:,2)==1&tbl_info(:,3)==1,[calcs_yy{i_calc_yy},'_1']} ;
    if yy_para==2 Scales =Scales./tbl_scale{tbl_info(:,2)==1&tbl_info(:,3)==1,[calcs_yy{i_calc_yy},'_1']}; end
    yy = Scales;
    mdl= fitlm(y,yy, 'linear','VarNames',[label_respondy,{label_respondyy}]);
    
    clear s;
    [~,~,~,~,s] = regress(yy,[ones(size(y)),y]);
    try stats.Coefficients = mdl.Coefficients; catch e, disp(['ERROR: ' e.message]); end  
    try stats.Rsquared = mdl.Rsquared; catch e, disp(['ERROR: ' e.message]); end
    stats.pValue = s(1,3);


    clear s;
    close all;
    figure('Unit','centimeter','Position',[5 5 10 10],'Visible','off'); hold on; 
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
    
    if stats.pValue>0.05
        sign = 'n.s.';
    elseif stats.pValue>0.01
        sign = '{\color{red}*}';
    elseif stats.pValue>0.001
        sign = '{\color{red}**}';
    else
        sign = '{\color{red}***}';
    end     
    s = sprintf('R^2=%.3f, adjusted R^2=%.3f, {\\itp}=%.4f (%s)', ...
     stats.Rsquared.Ordinary,stats.Rsquared.Adjusted,stats.pValue,sign);
    s = [upper(s(1)), s(2:end)];

    %s = strjoin(s,'\n');
    subtitle(s,'unit','normalized','Position',[0 1.05],'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8);
    set(gca,'FontSize',8,'LineWidth',1,'Unit','centimeter','Position',[1.5 4.5 7 4]);
    fn_fig = [label_respondy label_respondyy];
    fn_fig = replace(fn_fig,'\Delta','diff');
    

    if stats.pValue>0.05
        dn_save = [dn_results '\Preallnosign'];
    else
        dn_save = [dn_results '\Preallsign'];
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
Table=cell2table(tableall,'VariableNames',ColNames);

tablefile = '\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\CorrBold_individual_aseg\onePDF\presign.xls';
writetable(Table,tablefile,'Sheet','PreallFC_scale');
%% sperate stim and sham prebold-scale ： predict
clear mdl;
colorline = lines(length(calcs_yy));
ColNames = {'Stim','FC_para','FC','Sc_para','D-Scale','R2','Adjusted R2','p','Slope'};
i_stim=2;
tableall={};
k=1;
% Determination of ROI
%if ~exist(fullfile(fileparts(dn_results),'stat','MLR.txt'),'file') mkdir(fullfile(fileparts(dn_results),'stat','MLR.txt'));end
%dn_results = '\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\corrBOLD\NAc_EF_abs\'
diary(fullfile(fileparts(dn_results),'stat','Pre_individual_FC_aseg.txt'));
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
    y = tbl_outcomes{tbl_info(:,2)==1&tbl_info(:,3)==1,char(Y_para(y_para))};
    y = y(~(isnan(y)|isinf(y)));
    if isempty(y), continue; end

   % scale
    Scales = tbl_scale{tbl_info(:,2)==1&tbl_info(:,3)==1,[calcs_yy{i_calc_yy},'_2']}-tbl_scale{tbl_info(:,2)==1&tbl_info(:,3)==1,[calcs_yy{i_calc_yy},'_1']} ;
    if yy_para==2 Scales =Scales./tbl_scale{tbl_info(:,2)==1&tbl_info(:,3)==1,[calcs_yy{i_calc_yy},'_1']}; end
    yy = Scales;
    mdl{1} = fitlm(y(StimSubjs),yy(StimSubjs), 'linear','VarNames',[label_respondy,{label_respondyy}]);
    mdl{2} =fitlm(y(ShamSubjs),yy(ShamSubjs),'linear','VarNames',[label_respondy,{label_respondyy}]);
 
    clear s;
    [~,~,~,~,s(1,:)] = regress(yy(StimSubjs),[ones(size(y(StimSubjs))),y(StimSubjs)]);
    [~,~,~,~,s(2,:)] = regress(yy(ShamSubjs),[ones(size(y(ShamSubjs))),y(ShamSubjs)]);
    try stats(1).Coefficients = mdl{1}.Coefficients; catch e, disp(['ERROR: ' e.message]); end
    try stats(2).Coefficients = mdl{2}.Coefficients; catch e, disp(['ERROR: ' e.message]); end
    try stats(1).Rsquared = mdl{1}.Rsquared; catch e, disp(['ERROR: ' e.message]); end
    try stats(2).Rsquared = mdl{2}.Rsquared; catch e, disp(['ERROR: ' e.message]); end
    stats(1).pValue = s(1,3);
    stats(2).pValue = s(2,3);
    clear s;
    close all;
    figure('Unit','centimeter','Position',[5 5 10 10],'Visible','off'); hold on; 
    x_fit = 0.9*min(y):0.005:1.1*max(y);
    ylabel( label_respondyy);
    xlabel( label_respondy );
    for i_stim = length(stims):-1:1
        if i_stim == 1, c = colorline(i_calc_yy,:); pp = StimSubjs;
        else c = [0 0 0]; pp = ShamSubjs ;end
        x_stim = y(pp);
        y_stim = yy(pp);
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
     xlim([min(x_stim)-0.02*abs(min(x_stim)) max(x_stim)+0.02*abs(max(x_stim))]);


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
    fn_fig = [label_respondy label_respondyy];
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

tablefile = '\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\CorrBold_individual_aseg\onePDF\presign.xls';
writetable(Table,tablefile,'Sheet','PreFC_scale');
%% sperate stim and sham prebold+EF - scale ： predict（save sign to excle）
% Determination of ROI
%if ~exist(fullfile(fileparts(dn_results),'stat','MLR.txt'),'file') mkdir(fullfile(fileparts(dn_results),'stat','MLR.txt'));end
clear mdl;
colorline = lines(length(calcs_yy));
ColNames = {'Stim','EF_para','EF','FC_para','FC','Sc_para','D-Scale','R2','Adjusted R2','p','Slope'};
i_stim=2;
tableall={};
k=1;diary(fullfile(fileparts(dn_results),'stat','PreEFFC_scale_individual_aseg.txt'));
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
    y = tbl_outcomes{tbl_info(:,2)==1&tbl_info(:,3)==1,char(Y_para(y_para))};
    y = y(~(isnan(y)|isinf(y)));
    FC = y;
    if isempty(y), continue; end
    % EF
    EFvalues = tbl_EF{tbl_info(:,2)==1&tbl_info(:,3)==1,char(calcs_x(i_calc_x))} ;
    x = EFvalues;
    EF =x ;

   % scale
    Scales = tbl_scale{tbl_info(:,2)==1&tbl_info(:,3)==1,[calcs_yy{i_calc_yy},'_2']}-tbl_scale{tbl_info(:,2)==1&tbl_info(:,3)==1,[calcs_yy{i_calc_yy},'_1']} ;
    if yy_para==2 Scales =Scales./tbl_scale{tbl_info(:,2)==1&tbl_info(:,3)==1,[calcs_yy{i_calc_yy},'_1']}; end
    yy = Scales;
    % X
    load carsmall
    X = [EF,FC];
 
    %Table=table(EF,FC,Scales,'VariableNames',ColNames);

    % liner regression
    mdl{1} = fitlm(X(StimSubjs,:),yy(StimSubjs), 'linear','VarNames',[label_respondx,label_respondy,{label_respondyy}]);
    mdl{2} = fitlm(X(ShamSubjs,:),yy(ShamSubjs),'linear','VarNames',[label_respondx,label_respondy,{label_respondyy}]);
    % 
    clear s;
    X2 = [ones(size(EF)) EF FC EF.*FC];
    [b,BINT,R,RINT,STATS] = regress(yy(StimSubjs),X2(StimSubjs,:));
    [~,~,~,~,s(1,:)] = regress(yy(StimSubjs),X2(StimSubjs,:));
    [~,~,~,~,s(2,:)] = regress(yy(ShamSubjs),X2(ShamSubjs,:));
    try stats(1).Coefficients = mdl{1}.Coefficients; catch e, disp(['ERROR: ' e.message]); end
    try stats(2).Coefficients = mdl{2}.Coefficients; catch e, disp(['ERROR: ' e.message]); end
    try stats(1).Rsquared = mdl{1}.Rsquared; catch e, disp(['ERROR: ' e.message]); end
    try stats(2).Rsquared = mdl{2}.Rsquared; catch e, disp(['ERROR: ' e.message]); end
    stats(1).pValue = s(1,3);
    stats(2).pValue = s(2,3);
    clear s;
    close all;

   %figure('Unit','centimeter','Position',[5 5 10 10],'Visible','off'); hold on; 
 
    for i_stim = length(stims):-1:1
        if i_stim == 1, c = colorline(i_calc_yy,:);  pp = StimSubjs;
        else c = [0 0 0]; pp = ShamSubjs ;end
        x1=x(pp);
        x2 = y(pp);
        y_stim = yy(pp);
        scatter3(x1,x2,y_stim ,'filled','MarkerFaceColor',c);
        hold on  %在刚刚那副散点图上接着画
        if i_stim == 1
            x1fit = 0.9*min(x1):0.005:1.1*max(x1);;   %设置x1的数据间隔
            x2fit = 0.9*min(x2):0.005:1.1*max(x2);     %设置x2的数据间隔
            [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);    %生成一个二维网格平面，也可以说生成X1FIT,X2FIT的坐标
            hold on;
            YFIT=b(1)+b(2)*X1FIT+b(3)*X2FIT+b(4)*X1FIT.*X2FIT;    %代入已经求得的参数，拟合函数式
           mesh(X1FIT,X2FIT,YFIT);  
        end
        %view(10,10)  %改变角度观看已存在的三维图，第一个10表示方位角，第二个表示俯视角。
   

    end
    
    ylabel( label_respondy);
    xlabel( label_respondx );
    zlabel(label_respondyy);
    legend('Active','Sham','Location','southwest','FontSize',6);
     xlim([min(x)-0.02*abs(min(x)) max(x)+0.02*abs(max(x))]);


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
    fn_fig = [label_respondy label_respondyy];
    fn_fig = replace(fn_fig,'\Delta','diff');
    
    if stats(1).pValue>0.05
        dn_save = [dn_results '\PreEFFCnosign'];
    else
        dn_save = [dn_results '\PreEFFCsign'];
        i_stim=1;
        tableall(k,:) ={i_stim,X_para{x_para},label_respondx,Y_para{y_para},label_respondy,YY_para{yy_para}, label_respondyy,stats(1).Rsquared.Ordinary,...
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
end
end

diary off;
Table=cell2table(tableall,'VariableNames',ColNames);

tablefile = '\\172.16.4.19\eeg\Projects\tDCS\202111-SMHC_OCD\CorrBold_individual_aseg\onePDF\presign.xls';
writetable(Table,tablefile,'Sheet','PreEFFC_scale');
