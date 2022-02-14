%%
imgDir=dir('*.*');
sheet=2;
[ Num, Txt ] = xlsread( 'D:\jingzongfmri\stimTFnew.xlsx',sheet );
filename = Txt(2:46,[3,4]);
subjnum = length(filename);
%% get raw data from NAS
% 
DataP = '\\172.16.4.19\eeg\Projects\tDCS\202201-SMHC_cathode_FP2\data\fMRI\raw\post\WOC7-';
desPath = 'D:\jingzongfmri\data\fmri_indi_v2_new\';

for i=1:subjnum% 16
    desPathF = [desPath,'FunRaw\',char(filename(i,1)),'\'];
    desPathS = [desPath,'T1Raw\',char(filename(i,1)),'\'];
    if ~exist(desPathF,'file') mkdir(desPathF);end
    if ~exist(desPathS,'file') mkdir(desPathS);end
    dataFileF = [DataP,char(filename(i,2)),'*\*BOLD*'];
    dataFileS = [DataP,char(filename(i,2)),'*\*MPRAGE*'];
    Ff = dir(dataFileF);
    if length(Ff)==0 continue;end
    Sf = dir(dataFileS);
    fileF=[Ff.folder,'\',Ff.name];
   % fileS=[Sf(1).folder,'\',Sf(1).name];
    if exist(fileF,'file')
    copyfile( fileF,desPathF);
   % copyfile( fileS,desPathS);
    end    
end
%% get reorient from new data 

load('D:\jingzongfmri\data\fmri_indi_v1_new\get_reorient_from_newdata.mat')
Cfg.WorkingDir= 'D:\jingzongfmri\data\fmri_indi_v2_new';
Cfg.DataProcessDir='D:\jingzongfmri\data\fmri_indi_v2_new';
Cfg.SubjectID={'Sub_0041','Sub_0042','Sub_0043','Sub_0044','Sub_0045'}';
DPARSFA_run(Cfg);
%% get funImgARCFS 
load('D:\jingzongfmri\data\fmri_visit1_individual\Indi_1st_step_FunRaw.mat');
Cfg.WorkingDir= 'D:\jingzongfmri\data\fmri_indi_v1_new';
Cfg.DataProcessDir='D:\jingzongfmri\data\fmri_indi_v1_new';
Cfg.SubjectID={'Sub_0041','Sub_0042','Sub_0043','Sub_0044','Sub_0045'}';
DPARSFA_run(Cfg);
%% bet fun FunImgAR
% load('D:\jingzongfmri\data\fmri_visit1_individual\bet_FunImgAR.mat');
% Cfg.WorkingDir = 'D:\jingzongfmri\data\fmri_visit1_individual';
% Cfg.DataProcessDir = 'D:\jingzongfmri\data\fmri_visit1_individual';
% DPARSF_run(Cfg);
%% copy T1 indi  and  indi aspec to my folder
%DataP = '\\172.16.4.19\eeg\Projects\tDCS\202110-SMHC_OCD\recon\WOC7-';
%desPath = 'D:\jingzongfmri\data\fmri_visit1_individual\';
DataP = '\\172.16.4.19\eeg\Projects\tDCS\202201-SMHC_cathode_FP2\data\EF_simu\recon\WOC7-';
desPath = 'D:\jingzongfmri\data\fmri_indi_v1_new\';
%% copy T1 indi  and  indi aspec to my folder
for i=41:subjnum% 16
    desPathT = [desPath,'T1Indi\',char(filename(i,1)),'\'];
    desPathL = [desPath,'AspecIndi\',char(filename(i,1)),'\'];
    if ~exist(desPathT,'file') mkdir(desPathT);end
    if ~exist(desPathL,'file') mkdir(desPathL);end
    dataFileT = [DataP,char(filename(i,2)),'\outputs\brainmask.nii'];
    dataFileL = [DataP,char(filename(i,2)),'\outputs\label_area.nii'];
    if exist(dataFileT,'file')
        copyfile(dataFileT,desPathT);
        copyfile(dataFileL,desPathL);
    else 
        continue; %break;
    end
end
%% downsample 1mm to 3mm (spm)
for i=1:subjnum% 16
    desPathT = [desPath,'T1Indi\',char(filename(i,1)),'\brainmask.nii'];
    desPathL = [desPath,'AspecIndi\',char(filename(i,1)),'\label_area.nii'];
    if exist(desPathT,'file')
        lpz_resampleMRI(desPathT,[3 3 3])
        lpz_resampleMRI(desPathL,[3 3 3])
    else 
        continue; %break;
    end
end

%% copy funbet  to my folder
%BetP = 'D:\jingzongfmri\data\fmri_visit1_individual\RealignParameter\';
%desP = 'D:\jingzongfmri\data\fmri_visit1_individual\Betindicoreg\';
BetP = 'D:\jingzongfmri\data\fmri_indi_v2_new\RealignParameter\';
desP = 'D:\jingzongfmri\data\fmri_indi_v2_new\\Betindicoreg\';
%%
for i=41:subjnum% 16
    desPath = [desP,char(filename(i,1)),'\'];
    dataFile = [BetP,char(filename(i,1)),'\Bet_*.nii'];
    filed = dir(dataFile);
    if length(filed)
        if ~exist(desPath,'file') mkdir(desPath);end
        dfile = [filed.folder,'\',filed.name];
        copyfile(dfile,desPath);
    else 
        continue; %break;
    end
end
% %% Affine matrix applied to each time FunImg
% %% get affine matrix for each subj
% AfPath = 'D:\jingzongfmri\data\fmri_visit1_individual\T1Indi\';
% 
% for i=1:subjnum% 16
%     desPath = [desP,char(filename(i,1)),'\'];
%     if ~exist(desPath,'file') mkdir(desPath);end
%     dataFile = [AfPath,char(filename(i,1)),'\regf2a0GenericAffine.mat'];
%     Afmat = load(dataFile);
%     Afmat = Afmat.AffineTransform_double_3_3;
%     Afmat = reshape(Afmat,3,4);
%     Afmat(4,:)=[0,0,0,1]
%     filed = dir(dataFile);
%     dfile = [filed.folder,'\',filed.name];
%     if exist(dfile,'file')
%         copyfile(dfile,desPath);
%     else 
%         continue; %break;
%     end
% end
% 
% %% Bold process for 4 EF regions
% ROIsIntst = {'2012','2014','2027','2028'};
% %CorrROIsAllSubj=zeros(subjnum,2);


%% move FunImgARCFSold rsfMRI_reg to FunImgARCFS
%DataP = 'D:\jingzongfmri\data\fmri_visit2_individual\FunImgARCFSold\';
%desP = 'D:\jingzongfmri\data\fmri_visit2_individual\FunImgARCFS\';
DataP = 'D:\jingzongfmri\data\fmri_indi_v1_new\FunImgARCFSold\';
desP = 'D:\jingzongfmri\data\fmri_indi_v1_new\FunImgARCFS\';
%%
for i=41:subjnum% 16
    desPath = [desP,char(filename(i,1)),'\'];
   % Dataf = [DataP,char(filename(i,1)),'\*rsfMRI_reg.nii'];
    Dataf = [DataP,char(filename(i,1)),'\*rsfMRI_reg.nii'];
    filed = dir(Dataf);
    if length(filed)
        if ~exist(desPath,'file') mkdir(desPath);end
        dfile = [filed.folder,'\',filed.name];
        copyfile(dfile,desPath);
    else 
        continue; %break;
    end
end
%%subjnum = 39;
% %% copy T1Indi(NewT1) to T1Img
% DataP = '\\172.16.4.19\eeg\Projects\tDCS\202110-SMHC_OCD\recon\WOC7-';
% desP = 'D:\jingzongfmri\data\fmri_visit1_individual\T1Img\';
% %% copy T1 indi  and  indi aspec to my folder
% for i=1:subjnum% 16
%     desPath = [desP,char(filename(i,1)),'\'];
%     Dataf = [DataP,char(filename(i,2)),'\outputs\T1.nii'];
%     filed = dir(Dataf);
%     if length(filed)
%         if ~exist(desPath,'file') mkdir(desPath);end
%         dfile = [filed.folder,'\',filed.name];
%         desfile=[desPath,'T1',char(filename(i,1)),'.nii'];
%         copyfile(dfile,desfile);
%         %copyfile(dfile,desPath);
%     else 
%         continue; %break;
%     end
% end
% %% copy T1Indi(NewT1) to T1ImgCoreg
% DataP = '\\172.16.4.19\eeg\Projects\tDCS\202110-SMHC_OCD\recon\WOC7-';
% desP = 'D:\jingzongfmri\data\fmri_visit1_individual\T1ImgCoreg\';
% %% copy T1 indi  and  indi aspec to my folder
% for i=1:subjnum% 16
%     desPath = [desP,char(filename(i,1)),'\'];
%     Dataf = [DataP,char(filename(i,2)),'\outputs\T1.nii'];
%     filed = dir(Dataf);
%     if length(filed)
%         if ~exist(desPath,'file') mkdir(desPath);end
%         dfile = [filed.folder,'\',filed.name];
%         desfile=[desPath,'T1CoregFun',char(filename(i,1)),'.nii'];
%         copyfile(dfile,desfile);
%         %copyfile(dfile,desPath);
%     else 
%         continue; %break;
%     end
% end
% %% copy T1Indi(NewT1) to T1ImgBet
% DataP = 'D:\jingzongfmri\data\fmri_visit1_individual\T1Indi\';
% desP = 'D:\jingzongfmri\data\fmri_visit1_individual\T1ImgBet\';
% %%
% for i=1:subjnum% 16
%     desPath = [desP,char(filename(i,1)),'\'];
%     Dataf = [DataP,char(filename(i,1)),'\brainmask.nii'];
%     filed = dir(Dataf);
%     if length(filed)
%         if ~exist(desPath,'file') mkdir(desPath);end
%         dfile = [filed.folder,'\',filed.name];
%         desfile=[desPath,'Bet_T1',char(filename(i,1)),'.nii'];
%         copyfile(dfile,desfile);
%         %copyfile(dfile,desPath);
%     else 
%         continue; %break;
%     end
% end
% %% copy regWarped.nii.gz(NewFunbet) to RealignParameter\sub\Bet_meanaSub_001.nii
% DataP = 'D:\jingzongfmri\data\fmri_visit1_individual\FunImgARold\';
% desP = 'D:\jingzongfmri\data\fmri_visit1_individual\RealignParameter\';
% %%
% for i=1:subjnum% 16
%     desPath = [desP,char(filename(i,1)),'\'];
%     Dataf = [DataP,char(filename(i,1)),'\regWarped.nii.gz'];
%     filed = dir(Dataf);
%     if length(filed)
%         if ~exist(desPath,'file') mkdir(desPath);end
%         dfile = [filed.folder,'\',filed.name];
%         copyfile(dfile,desPath);
%         fold=dir([desPath,'Bet_mean*']);
%         fnew='regWarped.nii.gz';
%         if length(fold) 
%             cd(fold.folder)
%             new= ['old',fold.name];
%             old=fold.name ;
%             eval(['!rename' 32 old 32 new]); 
%             fold='regWarped.nii.gz';
%             fnew=['Bet_mean',char(filename(i,1)),'.nii']
%             eval(['!rename' , 32 fold 32 fnew]); 
%         end
%         
%     else 
%         continue; %break;
%     end
% end
%% get ROIs and FC from FunImgARCSF
% load Cfg
load('D:\jingzongfmri\data\fmri_indi_v1_new\getFCROIsIndi_2nd_step_FunImgARCFS_v1.mat')
%load('D:\jingzongfmri\data\fmri_visit1_individual\getFCROIsIndi_2nd_step_FunImgARCFS_v2.mat')
%load('D:\jingzongfmri\data\fmri_visit1_individual\Indi_2nd_step_FunImgARCFS.mat')
%load('D:\jingzongfmri\data\fmri_visit1_individual\T1newSeg_Cov_F008_S4_Indi.mat')
%%
%Cfg.StartingDirName = 'FunImgARCF01_08';
%Cfg.WorkingDir = 'D:\jingzongfmri\data\fmri_visit2_individual';
%Cfg.DataProcessDir = 'D:\jingzongfmri\data\fmri_visit2_individual';

ROIpath1='D:\jingzongfmri\data\fmri_indi_v1_new\AspecIndi\';
ROIpath2='\label_area.nii';
Maskpath1='D:\jingzongfmri\data\fmri_indi_v1_new\T1Indi\';
Maskpath2='\RS_brainmask.nii';
Cfg.MaskFile='';
Cfg.StartingDirName='FunImgARCFS';
Cfg.WorkingDir= 'D:\jingzongfmri\data\fmri_indi_v2_new';
Cfg.DataProcessDir='D:\jingzongfmri\data\fmri_indi_v2_new';
Cfg.SubjectID={'Sub_0041','Sub_0042','Sub_0043','Sub_0044','Sub_0045'}';
mySubjectID  = Cfg.SubjectID;
Cfg.SubjectID = [];
%%
for i =41:subjnum
    [bool_name,~] = ismember(filename(i,1), mySubjectID);
    if bool_name == 1
        Cfg.SubjectID = filename(i,1);
        Cfg.MaskFile=[Maskpath1,char(filename(i,1)),Maskpath2];

        Cfg.CalFC.ROIDef = {[ROIpath1,char(filename(i,1)),ROIpath2]};
        DPARSFA_run(Cfg);
    else
        continue;
    end 
end
