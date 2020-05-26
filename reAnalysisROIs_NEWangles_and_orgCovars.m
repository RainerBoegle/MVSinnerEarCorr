% This script contains all steps for the reAnalysis of the ROI data used in the
% JOON submission about inner ear parameter correlations with RS-fMRI activity.
% (reAnalysis because we made many changes and extensions in review round 1.)
%
% In short we use each of the angles ,e.g., left & right hSCC (horizontal semi-circular canals)
% and left & right Reid's plane (as measured from MR images), with the original covariates.
% ALSO, we will create the AVERAGE ANGLE angle from both sides.
% SPOILER WARNING ;) 
%     In the end we use pMVS and pReidsAv for reporting results.
%     (pMVS is the sinus-transformed angle of the average left and right hSCC angles.)
%     (pReidsAv is the sinus-transformed angle of the average left and right Reid's plane angles.)
%     In this script you can inspect all angles from left and right, 
%     AND of course the average angles that were reported in the manuscript.
%
% Specifically, we use each of the angles with the other covariates SEPARATELY,
% and compare the results in terms of the magnitude of R-squared (coefficient of
% determination) and contributions to R-squared (via stepwise regression) and
% the standard errors in the betas and the (normalized) betas of course.
% (Additionally AIC and BIC (Akaike Information Criterion and Bayesian
% Information Criterion) are used to compare between models.)
%
% The R-squared shows explanatory power for the given data, 
% standard error shows how well the estimated model can predict the data, or how
% uncertain a prediction would extected to be, e.g. in 95% of the cases is 2*se.
%
% All of this allows us to distinguish how well the TRANSFORMED angles are suited for
% "regressing out" the effects of MVS from RS-fMRI group data.
%
% IN THE END,
%    we want to know how good the pMVS (sinus-transformed hSCC angle), 
%    sinus-transformed Reids Plane for the left and right sides and average are 
%    for "regressing out" the effects of MVS from RS-fMRI group data.
%
%   This should allow researchers to make decisions about using CISS imaging or
%   just the Reid's plane angle for regressing out MVS effects on the group level.
%   SPOILER WARNING ;)
%     We suggest that CISS imaging should be used to extract the inner ear angle,
%     whenever available. Using Reid's plane is a "better than nothing" solution 
%     for older studies that do not have CISS imaging available.
%     There are some differences between left and right, but so far we deem these
%     not really noteworthy and only focus on the average angles.
%     (We have of course looked at the results from left and right, see below.)
%
%
%NB: the angles and parameters follow this naming convention
%     p*"Angle" is the sinus(Angle) transformed Angle of either
%        the Reid's Plane "pReids*"
%        the inner ear horizontal SCC "pMVS*"
%
%
%NB2: the other variables are Age, dt-RSfMRI, Sex and Handedness.
%     (dt-RSfMRI is the time spent in the MRI magnetic field at the start of RS-fMRI,
%     measured by the time difference of the start of the first sequence to the 
%     start of the RS-fMRI measurement)
%     For Age and dt-RSfMRI we also use the rank tranformation of the variables,
%     i.e. use tiedrank(X) on the variables as additional regressors.
%     All these "continuous" variables (everything besides Sex and Handedness)
%     are zscored before regression. 
%     The data is also zscored such that we get normalized betas from the regression.
%     Normalized betas show how many standard deviations of change in the data
%     are associated with a 1 standard deviation change in the regressors.
%     (E.G.: A beta of 0.75 for pMVS means that a 1 stdev change in pMVS changes the
%     RS-fMRI data by 0.75 stdev.)
%
%
%
%Author: Rainer.Boegle@googlemail.com
%V2.1: 26.05.2020
%Comment: V2.1:(26.05.2020): checked all comments and explanations, tried all plot modes. V2.0(20.05.2020): multiple bug fixes and adjustments for display of results and generation of tables. V1.0(13.05.2020): initial implementation, based on previous tests.


%% settings
plotLevel     = 'minimal'; %'minimal'; %'sparse'; %'all'; % %NB:'mininal'==all the plots for the manuscript + some extra. | 'sparse'==some more plots than in the manuscript. | 'all'==all the plots even just some diagnostic plotting of the covariates and data.
plotModeBetas = '2SE'; %plot 2*SE, i.e., 1*SE above and 1*SE below; 'SE'; %plot SE instead of CI, but just 1*SE above. And finally 'CI' for plotting the confidence intervals.
asOneFigure   = true;  %true==combine multiple data & trendlines figures by using subplots; choose false to see more details.
drop2Oldest   = false; %true; %SPOILER-WARNING: This does NOT change much. There are two subjects that are much older than all others, these can be removed, but that does not change much. Probably because we are regression out Age and tiedrank(Age)

%% data dir & path
dataDir     = pwd; %assume the data is in the same directory as the script file
fNameDataOI = 'dataReAnalysis.mat'; %the file name of the mat-file with the data struct.
pathDataOI  = [dataDir,filesep,fNameDataOI];

%% load the data
disp('loading the data of interest structure...');
dataReAnalysis = loadData(pathDataOI,drop2Oldest);

%% get basic infos about variables
disp('extracting basic infos & variables...');
roiNames = dataReAnalysis.roiNames;
nROIs = length(roiNames);

dataROIs = nanzscore(dataReAnalysis.dataROIs); %NB: we want normalized betas --> zscore data and regressors
nSubjs = size(dataROIs,1);

regressorsTable = dataReAnalysis.regressorsTable;

%% checking covariates & correlations
disp('checking covariates and (basic) correlations...');
[figsOutlier,figsCorr,corrStruct] = checkAllVariablesOI(dataROIs,roiNames,regressorsTable,plotLevel);

%% generate setup for Design Matrices for each transformed angle (l & r & av for pMVS & pReids)
baseDesMatSettings = {'Age_ALL','dt_RSfMRI_ALL','Sex','Handedness'}; 
baseDesMatString = curly(join(baseDesMatSettings,','),1);

%% run regression for each
disp(['running regression "Angle" & baseDesMat=="',baseDesMatString,'"...']);
disp('  Angle==pMVSL...');
fitResults_pMVSL = fitAllROIs(dataROIs,roiNames,regressorsTable,prepend('pMVSL',baseDesMatSettings));
disp('  Angle==pMVSR...');
fitResults_pMVSR = fitAllROIs(dataROIs,roiNames,regressorsTable,prepend('pMVSR',baseDesMatSettings));
disp('  Angle==pMVS...');
fitResults_pMVS  = fitAllROIs(dataROIs,roiNames,regressorsTable,prepend('pMVS', baseDesMatSettings));

disp('  Angle==pReidsL...');
fitResults_pReidsL = fitAllROIs(dataROIs,roiNames,regressorsTable,prepend('pReidsL', baseDesMatSettings));
disp('  Angle==pReidsR...');
fitResults_pReidsR = fitAllROIs(dataROIs,roiNames,regressorsTable,prepend('pReidsR', baseDesMatSettings));
disp('  Angle==pReidsAv...');
fitResults_pReidsAv= fitAllROIs(dataROIs,roiNames,regressorsTable,prepend('pReidsAv',baseDesMatSettings));


%% compare models using AIC and BIC
%% NB: lower is better
disp('checking if AIC or BIC of pMVS and pReidsAv are different...');
idx_pAngle = 1;
idx_other= 2:size(fitResults_pMVSL.desMat,2);
pVals = zeros(nROIs,1);
tVals = zeros(nROIs,1);
dfe   = zeros(nROIs,1);
disp('  The difference in AIC or BIC for pMVS VS pReidsAv,...');
charROInames = char(roiNames);
for indROI = 1:nROIs
    %AIC
    currAIC_pMVS = fitResults_pMVS.AIC( indROI);
    currAIC_pReidsAv= fitResults_pReidsAv.AIC(indROI);
    if(currAIC_pMVS < currAIC_pReidsAv)
        infoStrAIC = ['AIC is BETTER (lower) for "pMVS" than "pReidsAv". (',num2str(currAIC_pMVS),'<',num2str(currAIC_pReidsAv),')'];
    else
        infoStrAIC = ['AIC is worse (HIGHER) for "pMVS" than "pReidsAv". (',num2str(currAIC_pMVS),'>',num2str(currAIC_pReidsAv),')'];
    end
    %BIC
    currBIC_pMVS  = fitResults_pMVS.BIC(indROI);
    currBIC_pReidsAv = fitResults_pReidsAv.BIC(indROI);
    if(currBIC_pMVS < currBIC_pReidsAv)
        infoStrBIC = ['BIC is BETTER (lower) for "pMVS" than "pReidsAv". (',num2str(currBIC_pMVS),'<',num2str(currBIC_pReidsAv),')'];
    else
        infoStrBIC = ['BIC is worse (HIGHER) for "pMVS" than "pReidsAv". (',num2str(currBIC_pMVS),'>',num2str(currBIC_pReidsAv),')'];
    end
    
    disp(['    at the ROI "',charROInames(indROI,:),'" is:']);
    disp(['                ',infoStrAIC]);
    disp(['                ',infoStrBIC]);
end

%% show main results of regression for each setup 
% betas and StdErr for each regressor 
% (StdErr instead of CI, because StdErr shows how good the model is for predicting)
% In the code I write SE rather than StdErr.

disp('displaying robust regression results...');
switch(plotLevel)
    case 'all'
        disp('  Angle==pMVSL...');
        figsRegResults_pMVSL = plotFitResults(    fitResults_pMVSL,4000+420,'pMVSL',plotModeBetas,[],plotLevel);
        figsDataTrends_pMVSL = plotDataTrendlines(fitResults_pMVSL,idx_pAngle,idx_other,1000+230,asOneFigure);
        disp('  Angle==pMVSR...');
        figsRegResults_pMVSR = plotFitResults(    fitResults_pMVSR,3000+420,'pMVSR',plotModeBetas,[],plotLevel);
        figsDataTrends_pMVSR = plotDataTrendlines(fitResults_pMVSR,idx_pAngle,idx_other,2000+230,asOneFigure);
        disp('  Angle==pMVS...');
        figsRegResults_pMVS  = plotFitResults(    fitResults_pMVS, 1000+420,'pMVS', plotModeBetas,[],plotLevel);
        figsDataTrends_pMVS  = plotDataTrendlines(fitResults_pMVS, idx_pAngle,idx_other,3000+230,asOneFigure);
    case 'sparse'
        disp('  Angle==pMVS...');
        figsRegResults_pMVS  = plotFitResults(    fitResults_pMVS, 1000+420,'pMVS', plotModeBetas,[],plotLevel);
        figsDataTrends_pMVS  = plotDataTrendlines(fitResults_pMVS, idx_pAngle,idx_other,3000+230,asOneFigure);
    case 'minimal'
        disp('  Angle==pMVS...');
        figsRegResults_pMVS  = plotFitResults(    fitResults_pMVS, 1000+420,'pMVS', plotModeBetas,'pMVS',plotLevel);
        figsDataTrends_pMVS  = plotDataTrendlines(fitResults_pMVS, idx_pAngle,idx_other,3000+230,asOneFigure);
    otherwise
        error(['ERROR: unkown plot level "',plotLevel,'"!']);
end

switch(plotLevel)
    case 'all'
        disp('  Angle==pReidsL...');
        figsRegResults_pReidsL  = plotFitResults(    fitResults_pReidsL,5000+420,'pReidsL',plotModeBetas,[],plotLevel);
        figsDataTrends_pReidsL  = plotDataTrendlines(fitResults_pReidsL,idx_pAngle,idx_other,4000+230,asOneFigure);
        disp('  Angle==pReidsR...');
        figsRegResults_pReidsR  = plotFitResults(    fitResults_pReidsR,6000+420,'pReidsR',plotModeBetas,[],plotLevel);
        figsDataTrends_pReidsR  = plotDataTrendlines(fitResults_pReidsR,idx_pAngle,idx_other,5000+230,asOneFigure);
        disp('  Angle==pReidsAv...');
        figsRegResults_pReidsAv = plotFitResults(    fitResults_pReidsAv,2000+420,'pReidsAv',plotModeBetas,[],plotLevel);
        figsDataTrends_pReidsAv = plotDataTrendlines(fitResults_pReidsAv,idx_pAngle,idx_other,6000+230,asOneFigure);
    case 'sparse'
        disp('  Angle==pReidsAv...');
        figsRegResults_pReidsAv = plotFitResults(    fitResults_pReidsAv,2000+420,'pReidsAv',plotModeBetas,[],plotLevel);
        figsDataTrends_pReidsAv = plotDataTrendlines(fitResults_pReidsAv,idx_pAngle,idx_other,6000+230,asOneFigure);
    case 'minimal'
        disp('  Angle==pReidsAv...');
        figsRegResults_pReidsAv = plotFitResults(    fitResults_pReidsAv,2000+420,'pReidsAv',plotModeBetas,'pReidsAv',plotLevel);
        figsDataTrends_pReidsAv = plotDataTrendlines(fitResults_pReidsAv,idx_pAngle,idx_other,6000+230,asOneFigure);
    otherwise
        error(['ERROR: unkown plot level "',plotLevel,'"!']);
end
drawnow;   

%% prepare & run stepwise regression for each setup
disp('preparing & running stepwise regression for each setup...');
combineSexHandedness = false; %true; %--> allow sex & handedness as factors but if both are there, replace by Sex*Handedness interaction (avoids making the model rank deficient)
addConst = true; %false; %Without a constant term there can be models with negative R-squared as the whole model is non-negative in R-square but some of the reduced models of the initial steps especially step 1 with the first regressor can be problematic, therefore we add a constant.

disp('  Angle==pMVS...');
resultsStepwiseReg_pMVS     = prepAndRunStepwiseReg(prepend('pMVS',    baseDesMatSettings),dataROIs,roiNames,regressorsTable,combineSexHandedness,addConst);
disp('  Angle==pReidsAv...');
resultsStepwiseReg_pReidsAv = prepAndRunStepwiseReg(prepend('pReidsAv',baseDesMatSettings),dataROIs,roiNames,regressorsTable,combineSexHandedness,addConst);

if(strcmp(plotLevel,'all')) %stepwise regression takes a loooong time therefore I don't bother with it if I don't display it later.
disp('  Angle==pMVSL...');
resultsStepwiseReg_pMVSL    = prepAndRunStepwiseReg(prepend('pMVSL',  baseDesMatSettings),dataROIs,roiNames,regressorsTable,combineSexHandedness,addConst);
disp('  Angle==pMVSR...');
resultsStepwiseReg_pMVSR    = prepAndRunStepwiseReg(prepend('pMVSR',  baseDesMatSettings),dataROIs,roiNames,regressorsTable,combineSexHandedness,addConst);
disp('  Angle==pReidsL...');
resultsStepwiseReg_pReidsL  = prepAndRunStepwiseReg(prepend('pReidsL',baseDesMatSettings),dataROIs,roiNames,regressorsTable,combineSexHandedness,addConst);
disp('  Angle==pReidsR...');
resultsStepwiseReg_pReidsR  = prepAndRunStepwiseReg(prepend('pReidsR',baseDesMatSettings),dataROIs,roiNames,regressorsTable,combineSexHandedness,addConst);
end


%% display each stepwise regression results for each setup
disp('displaying each stepwise regression results for each setup...');
rSquareYLimits = [0 max([60; max([floor(max(squeeze(sum(resultsStepwiseReg_pMVS.addedRsquare)).*100)*1.1); floor(max(squeeze(sum(resultsStepwiseReg_pReidsAv.addedRsquare)).*100)*1.1)])])];

disp('  Angle==pMVS...');
figStepwiseReg_pMVS     = plotStepwiseRegResults(resultsStepwiseReg_pMVS,    421,rSquareYLimits);
disp('  Angle==pReidsAv...');
figStepwiseReg_pReidsAv = plotStepwiseRegResults(resultsStepwiseReg_pReidsAv,422,rSquareYLimits);

if(strcmp(plotLevel,'all'))
disp('  Angle==pMVSL...');
figStepwiseReg_pMVSL    = plotStepwiseRegResults(resultsStepwiseReg_pMVSL,   423,rSquareYLimits);
disp('  Angle==pMVSR...');
figStepwiseReg_pMVSR    = plotStepwiseRegResults(resultsStepwiseReg_pMVSR,   424,rSquareYLimits);
disp('  Angle==pReidsL...');
figStepwiseReg_pReidsL  = plotStepwiseRegResults(resultsStepwiseReg_pReidsL, 425,rSquareYLimits);
disp('  Angle==pReidsR...');
figStepwiseReg_pReidsR  = plotStepwiseRegResults(resultsStepwiseReg_pReidsR, 426,rSquareYLimits);
end


%% make tables (format results for table in the manuscript --> some more details in this output than in the manuscript as we had to make it smaller in the manuscript.)
%% NB (why I show SEs in review round 1, rather than CIs as in the original manuscript):
% I now rather want betas and stdErr of the model per region to be displayed and reported.
% This is better for comparing normalized beta as MODULATION STRENGTH across transformed angles
% (normalized betas allow to say how many standard deviations (Std) of change happen 
%  ,e.g., for beta==0.5 we have a change in 0.5*Std of the data for every 1 Std in the regressor.)
% and
% to compare stdErr size as PRECISION of prediction 
%  --> https://blog.minitab.com/blog/adventures-in-statistics-2/regression-analysis-how-to-interpret-s-the-standard-error-of-the-regression
%  also keep in mind --> https://blog.minitab.com/blog/adventures-in-statistics-2/how-high-should-r-squared-be-in-regression-analysis
%  and --> https://blog.minitab.com/blog/adventures-in-statistics-2/applied-regression-analysis-how-to-present-and-use-the-results-to-avoid-costly-mistakes-part-2
%  OR  --> https://blog.minitab.com/blog/adventures-in-statistics-2/when-should-i-use-confidence-intervals-prediction-intervals-and-tolerance-intervals
%  also consider this comment --> https://blog.minitab.com/blog/adventures-in-statistics-2/multiple-regession-analysis-use-adjusted-r-squared-and-predicted-r-squared-to-include-the-correct-number-of-variables
%
%From the above:
%"Check the Standard Error of the Regression (SE): 
%    R-squared gets all of the attention; however, it does not tell you how 
%    the data values compare to the predicted values. SE does just that! 
%    SE is measured in the units of the response variable and represents 
%    the standard distance data values fall from the regression line. 
%    Also, the fitted value plus/minus 2*standard error of the regression 
%    provides a quick approximation of a 95% prediction interval."#

disp('creating tables from regression results...');
outXLSfileName = ['reAnalysisResults_',datestr(now,'yyyy_mmm_dd_HHMM'),'.xlsx'];

%% one sheet for pValues 
legendStrRegs    = fitResults_pMVS.colNames_desMat(:);
legendStrRegs{1} = 'pMVS | pReidsAv';
pValJoinFun = @(X,Y) cellfun(@(pValStr1,pValStr2) [pValStr1,' | ',pValStr2],formatPVal( X ),formatPVal( Y ),'UniformOutput',false);
table_pValues = table(pValJoinFun(fitResults_pMVS.pVals(:,1),fitResults_pReidsAv.pVals(:,1)),...
                      pValJoinFun(fitResults_pMVS.pVals(:,2),fitResults_pReidsAv.pVals(:,2)),...
                      pValJoinFun(fitResults_pMVS.pVals(:,3),fitResults_pReidsAv.pVals(:,3)),...
                      pValJoinFun(fitResults_pMVS.pVals(:,4),fitResults_pReidsAv.pVals(:,4)),...
                      pValJoinFun(fitResults_pMVS.pVals(:,5),fitResults_pReidsAv.pVals(:,5)),...
                      pValJoinFun(fitResults_pMVS.pVals(:,6),fitResults_pReidsAv.pVals(:,6)),...
                      'VariableNames',roiNames,'RowNames',legendStrRegs);
writetable(table_pValues,outXLSfileName,'WriteRowNames',true,'Sheet','pValues');

%% one sheet for betas & CI or SE (default)
switch(plotModeBetas)
    case 'CI'
        %% betas 
        table_betas = table(fitResults_pMVS.betas(:,1),fitResults_pMVS.betas(:,2),fitResults_pMVS.betas(:,3),...
                            fitResults_pMVS.betas(:,4),fitResults_pMVS.betas(:,5),fitResults_pMVS.betas(:,6),...
                            'VariableNames',roiNames,'RowNames',fitResults_pMVS.colNames_desMat(:));
        writetable(table_betas,outXLSfileName,'WriteRowNames',true,'Sheet','betas - pMVS');
        
        table_betas = table(fitResults_pReidsAv.betas(:,1),fitResults_pReidsAv.betas(:,2),fitResults_pReidsAv.betas(:,3),...
                            fitResults_pReidsAv.betas(:,4),fitResults_pReidsAv.betas(:,5),fitResults_pReidsAv.betas(:,6),...
                            'VariableNames',roiNames,'RowNames',fitResults_pReidsAv.colNames_desMat(:));
        writetable(table_betas,outXLSfileName,'WriteRowNames',true,'Sheet','betas - pReidsAv');

        %% betaIs
        table_betaIs = table(fitResults_pMVS.betaIs(:,:,1)',fitResults_pMVS.betaIs(:,:,2)',fitResults_pMVS.betaIs(:,:,3)',...
                             fitResults_pMVS.betaIs(:,:,4)',fitResults_pMVS.betaIs(:,:,5)',fitResults_pMVS.betaIs(:,:,6)',...
                             'VariableNames',roiNames,'RowNames',fitResults_pMVS.colNames_desMat(:));
        writetable(table_betaIs,outXLSfileName,'WriteRowNames',true,'Sheet','beta CIs - pMVS');
        
        table_betaIs = table(fitResults_pReidsAv.betaIs(:,:,1)',fitResults_pReidsAv.betaIs(:,:,2)',fitResults_pReidsAv.betaIs(:,:,3)',...
                             fitResults_pReidsAv.betaIs(:,:,4)',fitResults_pReidsAv.betaIs(:,:,5)',fitResults_pReidsAv.betaIs(:,:,6)',...
                             'VariableNames',roiNames,'RowNames',fitResults_pReidsAv.colNames_desMat(:));
        writetable(table_betaIs,outXLSfileName,'WriteRowNames',true,'Sheet','beta CIs - pReidsAv');
        
    case {'SE';'2SE'}
        %% betas & SE together
%         table_BetasSEs = table(joinResults2Str(fitResults_pMVS.betas(:,1),fitResults_pMVS.stdErr(:,1)),joinResults2Str(fitResults_pMVS.betas(:,2),fitResults_pMVS.stdErr(:,2)),joinResults2Str(fitResults_pMVS.betas(:,3),fitResults_pMVS.stdErr(:,3)),...
%                                joinResults2Str(fitResults_pMVS.betas(:,4),fitResults_pMVS.stdErr(:,4)),joinResults2Str(fitResults_pMVS.betas(:,5),fitResults_pMVS.stdErr(:,5)),joinResults2Str(fitResults_pMVS.betas(:,6),fitResults_pMVS.stdErr(:,6)),...
%                                'VariableNames',roiNames,'RowNames',fitResults_pMVS.colNames_desMat(:));
%         writetable(table_BetasSEs,outXLSfileName,'WriteRowNames',true,'Sheet','betas SEs - pMVS');
%         
%         table_BetasSEs = table(joinResults2Str(fitResults_pReidsAv.betas(:,1),fitResults_pReidsAv.stdErr(:,1)),joinResults2Str(fitResults_pReidsAv.betas(:,2),fitResults_pReidsAv.stdErr(:,2)),joinResults2Str(fitResults_pReidsAv.betas(:,3),fitResults_pReidsAv.stdErr(:,3)),...
%                                joinResults2Str(fitResults_pReidsAv.betas(:,4),fitResults_pReidsAv.stdErr(:,4)),joinResults2Str(fitResults_pReidsAv.betas(:,5),fitResults_pReidsAv.stdErr(:,5)),joinResults2Str(fitResults_pReidsAv.betas(:,6),fitResults_pReidsAv.stdErr(:,6)),...
%                                'VariableNames',roiNames,'RowNames',fitResults_pReidsAv.colNames_desMat(:));
%         writetable(table_BetasSEs,outXLSfileName,'WriteRowNames',true,'Sheet','betas SEs - pReidsAv');
        
        
        legendStrRegs = regexprep(regexprep(fitResults_pMVS.colNames_desMat(:),'_','-'),'-ALL','');
        legendStrRegs{1} = 'pMVS | pReidsAv';
        joinStrAdj = @(Betas1,SEs1,Betas2,SEs2) cellfun(@(str1,str2) [str1,' | ',str2],joinResults2Str(Betas1,SEs1),joinResults2Str(Betas2,SEs2),'UniformOutput',false);
        table_BetasSEs = table(joinStrAdj(fitResults_pMVS.betas(:,1),fitResults_pMVS.stdErr(:,1),fitResults_pReidsAv.betas(:,1),fitResults_pReidsAv.stdErr(:,1)),...
                               joinStrAdj(fitResults_pMVS.betas(:,2),fitResults_pMVS.stdErr(:,2),fitResults_pReidsAv.betas(:,2),fitResults_pReidsAv.stdErr(:,2)),...
                               joinStrAdj(fitResults_pMVS.betas(:,3),fitResults_pMVS.stdErr(:,3),fitResults_pReidsAv.betas(:,3),fitResults_pReidsAv.stdErr(:,3)),...
                               joinStrAdj(fitResults_pMVS.betas(:,4),fitResults_pMVS.stdErr(:,4),fitResults_pReidsAv.betas(:,4),fitResults_pReidsAv.stdErr(:,4)),...
                               joinStrAdj(fitResults_pMVS.betas(:,5),fitResults_pMVS.stdErr(:,5),fitResults_pReidsAv.betas(:,5),fitResults_pReidsAv.stdErr(:,5)),...
                               joinStrAdj(fitResults_pMVS.betas(:,6),fitResults_pMVS.stdErr(:,6),fitResults_pReidsAv.betas(:,6),fitResults_pReidsAv.stdErr(:,6)),...
                               'VariableNames',roiNames,'RowNames',legendStrRegs);
        writetable(table_BetasSEs,outXLSfileName,'WriteRowNames',true,'Sheet','betas SEs');
        
    otherwise
        error(['ERROR: unkown plot mode "',plotModeBetas,'".']);
end

%% one sheet for added R-square
legendStrRegs = regexprep(regexprep(resultsStepwiseReg_pMVS.desMatSettings,'_','-'),'-ALL','');
legendStrRegs{1} = 'pMVS | pReidsAv';
joinResults2StrAdj = @(x,y) joinResults2Str(x,y,' | ',2,2,false);
addedRsquarePerc_pMVS     = round(resultsStepwiseReg_pMVS.addedRsquare.*10000)./100;
addedRsquarePerc_pReidsAv = round(resultsStepwiseReg_pReidsAv.addedRsquare.*10000)./100;
table_addedRsquare = table(joinResults2StrAdj(addedRsquarePerc_pMVS(:,1),addedRsquarePerc_pReidsAv(:,1)),joinResults2StrAdj(addedRsquarePerc_pMVS(:,2),addedRsquarePerc_pReidsAv(:,2)),joinResults2StrAdj(addedRsquarePerc_pMVS(:,3),addedRsquarePerc_pReidsAv(:,3)),...
                           joinResults2StrAdj(addedRsquarePerc_pMVS(:,4),addedRsquarePerc_pReidsAv(:,4)),joinResults2StrAdj(addedRsquarePerc_pMVS(:,5),addedRsquarePerc_pReidsAv(:,5)),joinResults2StrAdj(addedRsquarePerc_pMVS(:,6),addedRsquarePerc_pReidsAv(:,6)),...
                           'VariableNames',roiNames,'RowNames',legendStrRegs);
writetable(table_addedRsquare,outXLSfileName,'WriteRowNames',true,'Sheet','addedRsquare');

%% done
disp(' ');
disp('Done.');
disp(' ');

%% subfunction 
%% plotStepwiseResults
function figObj = plotStepwiseRegResults(resultsStepwiseReg,figNum,rSquareYLimits)
% plot results of stepwise regression

addedRsquare     = resultsStepwiseReg.addedRsquare;
fracAddedRsquare = resultsStepwiseReg.fracAddedRsquare;
roiNames         = resultsStepwiseReg.roiNames;
nROIs = length(roiNames);
legendStrRegs = regexprep(regexprep(resultsStepwiseReg.desMatSettings,'_','-'),'-ALL','');
if(~exist('rSquareYLimits','var')||isempty(rSquareYLimits))
    rSquareYLimits = [0 max([floor(max(squeeze(sum(addedRsquare)).*100)*1.1); 60])];
end

%% final plot: stacked bar graph with median Rsquared per ROI
figObj = figure(figNum); clf;
subplot(2,1,1);
bar(1:nROIs,addedRsquare'.*100,'stacked'); title('added Rsquare (%-explained variance added by each regressor for each ROI)');
xticks(1:nROIs); xticklabels(roiNames);
yticks(0:10:60);
legend(legendStrRegs);
%ylabel('added %-explained variance');
ylim(rSquareYLimits)

subplot(2,1,2);
bar(1:nROIs,fracAddedRsquare','stacked'); title('fraction of added Rsquare relative to total Rsquare');
xticks(1:nROIs); xticklabels(roiNames);
legend(legendStrRegs);
%ylabel('fraction of added %-explained variance relative to total');
ylim([0, 1.1]);

end

%% prepAndRunStepwiseReg
function resultsStepwiseReg = prepAndRunStepwiseReg(desMatSettings,dataROIs,roiNames,regressorsTable,combineSexHandedness,addConst)
% do the full steps of the stepwise regression appraoch, i.e., 
% create all permutations
% run stepwise regression with all permutations, adding each regressor one at a time.
% take all R-squared values and determine the differences --> the addedRsquare for each time a regressor was added.
% generate the median R-squared values and the fractions

nROIs = size(dataROIs,2);
nRegs = length(desMatSettings);

%% create all permutations
disp('   setting up permutations for stepwise (robust) regression...')
allPermsRegressors = perms(desMatSettings);
nSteps = size(allPermsRegressors,1);

%% run stepwise regression
disp('   running stepwise (robust) regression...');
allRsquare = runStepwiseRegression(allPermsRegressors,dataROIs,roiNames,regressorsTable,addConst,combineSexHandedness);

%% generate the addedRsquare
disp('   determining added Rsquare from stepwise regression results...');
[addedRsquare,allAddedRsquare] = determineAddedRsquare(allRsquare,allPermsRegressors,desMatSettings);
allFracAddedRsquare = zeros(size(allAddedRsquare));
for indROI = 1:nROIs
    for indStep = 1:nSteps
        currData = squeeze(allAddedRsquare(indStep,:,indROI));
        allFracAddedRsquare(indStep,:,indROI) = currData./sum(currData);
    end
end
medianFracAddedRsquare = squeeze(median(allFracAddedRsquare,1));
medianFracAddedRsquare = medianFracAddedRsquare./sum(medianFracAddedRsquare,1);
fracAddedRsquare = (addedRsquare./repmat(sum(addedRsquare),nRegs,1));

%% create output struct
resultsStepwiseReg = struct('desMatSettings',{desMatSettings},...
                            'roiNames',{roiNames},...
                            'allPermsRegressors',{allPermsRegressors},...
                            'allRsquare',allRsquare,...
                            'allAddedRsquare',allAddedRsquare,...
                            'addedRsquare',addedRsquare,...
                            'allFracAddedRsquare',allFracAddedRsquare,...
                            'medianFracAddedRsquare',medianFracAddedRsquare,...
                            'fracAddedRsquare',fracAddedRsquare);

end

%% determineAddedRsquare
function [addedRsquare,allAddedRsquare] = determineAddedRsquare(allRsquare,allPermsRegressors,baseStepwiseReg)
% determine the added Rsquare for each type of regressors

nSteps = size(allPermsRegressors,1);
nRegsTotal = size(allPermsRegressors,2);
nROIs = size(allRsquare,3);
allAddedRsquare = zeros(nSteps,nRegsTotal,nROIs);
for indStep = 1:nSteps
    for indRegMax = 1:nRegsTotal
        if(indRegMax==1)
            refRsquare = zeros(nROIs,1);
        else
            refRsquare = squeeze(allRsquare(indStep,indRegMax-1,:));
        end
        currRegName = allPermsRegressors{indStep,indRegMax};
        assignIndex = find(cellfun(@(str) strcmp(currRegName,str),baseStepwiseReg));
        if(length(assignIndex)~=1)
            error('assignIndex should only be one index not multiple or empty!!! Check what is going on! This should not happen!');
        end
        
        allAddedRsquare(indStep,assignIndex,:) = squeeze(allRsquare(indStep,indRegMax,:)) - refRsquare;
    end
end
addedRsquare = squeeze(median(allAddedRsquare,1));
end

%% runStepwiseRegression
function allRsquare = runStepwiseRegression(allPermsRegressors,dataROIs,roiNames,regressorsTable,addConst,combineSexHandedness)
% run the stepwise regression and get Rsquare values back 
% (later check values and create median Rsquare difference with different function)

nSteps = size(allPermsRegressors,1);
nRegsTotal = size(allPermsRegressors,2);
nROIs = size(dataROIs,2);
allRsquare = zeros(nSteps,nRegsTotal,nROIs);
reverseStr = ''; %init
for indStep = 1:nSteps
    %disp(['Step ',num2str(indStep),' of ',num2str(nSteps)]);
    msg = sprintf('    Step %d of %d', indStep, nSteps); 
    fprintf([reverseStr, msg]); %delete the last message and print current message.
    reverseStr = repmat(sprintf('\b'), 1, length(msg)); %make current message length of backspace for next deletion
    for indRegMax = 1:nRegsTotal
        %settings for design matrix
        if(addConst)
            %disp(['  Regs 1:',num2str(indRegMax),' + "const"...']);
            currDesMatSettings = allPermsRegressors(indStep,1:indRegMax);
            currDesMatSettings{end+1} = 'const';
        else
            %disp(['  Regs 1:',num2str(indRegMax),'...']);
            currDesMatSettings = allPermsRegressors(indStep,1:indRegMax);
        end
        if(combineSexHandedness)
            currDesMatSettings = replaceSexHandedness(currDesMatSettings);
        end
        %currDesMatSettings = replaceALLtypes(currDesMatSettings);
        
        fitResults = fitAllROIs(dataROIs,roiNames,regressorsTable,currDesMatSettings);
        currRsquare = fitResults.stats(1,:)';
        if(any(currRsquare<0))
            disp(['WARNING: current Rsquare is negative!!! (Rsq==',num2str(currRsquare(:)'),')']);
            keyboard; %I will think about the way to solve this later, e.g., always use const when using one regressor or always adding a const in this case......
        end
        allRsquare(indStep,indRegMax,:) = currRsquare;
    end
end
disp(' ');
end

%% replaceSexHandedness
function newDesMatSettings = replaceSexHandedness(desMatSettings)
% check if settings has Sex & Handedness, if yes, replace by Sex*Handedness


selSex = cellfun(@(str) strcmp(str,'Sex'),desMatSettings);
hasSex = any(selSex);
selHandedness = cellfun(@(str) strcmp(str,'Handedness'),desMatSettings);
hasHandedness = any(selHandedness);
if(hasSex&&hasHandedness)
    selSexHandedness = selSex|selHandedness;
    selOthers = ~selSexHandedness;
    if(any(selOthers))
        newDesMatSettings = desMatSettings(selOthers); %keep all other regressors
        newDesMatSettings{end+1} = 'Sex*Handedness';
    else %there are no other regressors besides Sex & Handedness
        newDesMatSettings = {'Sex*Handedness'}; 
    end
else
    newDesMatSettings = desMatSettings; %leave unchanged
end

end

%% joinBetaSE
function csBetasSEs = joinResults2Str(in1,in2,joinStr,precisionN1,precisionN2,addSign)
% join values as strings with extra markings 
if(~exist('joinStr','var')||isempty(joinStr))
    joinStr = '\pm';
end
if(~exist('precisionN1','var')||isempty(precisionN1))
    precisionN1 = 3;
end
if(~exist('precisionN2','var')||isempty(precisionN2))
    precisionN2 = 3;
end
if(~exist('addSign','var')||isempty(addSign))
    addSign = true;
end


csBetasSEs = cell(length(in1),1);
for ind = 1:length(in1)
    if(in1(ind)>0&&addSign)
        signStr = '+';
        betaStr = abs(in1(ind));
    elseif(in1(ind)<0&&addSign)
        signStr = '-';
        betaStr = abs(in1(ind));
    else
        signStr = '';
        betaStr = in1(ind);
    end
    csBetasSEs{ind} = [signStr,num2strN(betaStr,precisionN1),' ',joinStr,' ',num2strN(in2(ind),precisionN2)];
end
end

%% num2strN
function outStr = num2strN(num,nDigits)
% num2str but limit nDigits

if(mod(num,1)==0)
    outStr = num2str(num);
else
    outStr = num2str(num);
    nMax = length(outStr);
    idxDot = regexp(outStr,'\.');
    if(length(outStr)>idxDot)
        outStr = outStr([1:idxDot,(idxDot+1):min([nMax; (idxDot+nDigits)])]);
    end
end
    
end
%% formatPVal
function csPVals = formatPVal(pVals)
% format all pVals and put in cellstring

csPVals = cell(length(pVals),1);
for ind = 1:length(pVals)
    csPVals{ind} = formatPValStr(pVals(ind));
end
end

%% formatPValStr
function pValStr = formatPValStr(pVal)
% format the pValue string 

if(pVal<0.001)
    pValStr = num2str(pVal,'%3.2e');
elseif((pVal>=0.001)&&(pVal<=(0.05/6)))
    pValStr = num2str(pVal,'%5.4f');
else
    pValStr = num2str(pVal,'%4.3f');
end
if(str2num(pValStr(end))==0)
    pValStr(end) = [];
end

end

%% plotDataTrendlines
function figsDataTrends = plotDataTrendlines(fitResults,colOI,colsRem,figBaseNum,asOneFig)
% plot data & trendlines for ColOI --> first remove effects of colsRem

if(~exist('figBaseNum','var')||isempty(figBaseNum))
    figBaseNum = 230;
end

if(~exist('asOneFig','var')||isempty(asOneFig))
    asOneFig = false;
end

%% get data
dataROIs = fitResults.dataROIs;
roiLabels = regexprep(fitResults.roiNames,'_','\\_');
nROIs = size(dataROIs,2);
desMatRem = fitResults.desMat(:,colsRem);
betasRem  = fitResults.betas(   colsRem,:);
joinNameRem = join(fitResults.colNames_desMat(colsRem),';');
remStr = joinNameRem{1};

regColOI = squeeze(fitResults.desMat(:,colOI));
colNameOI= fitResults.colNames_desMat{ colOI};
betaOI   = squeeze(fitResults.betas(   colOI,:));
trendX   = linspace(min(regColOI),max(regColOI),100);

%% plot
if(asOneFig)
    figsDataTrends = cell(nROIs+1,1);
else
    figsDataTrends = cell(nROIs,1);
end
for indROI = 1:nROIs
    if(asOneFig)
        if(indROI==1)
            figsDataTrends{1} = figure(figBaseNum); clf;
        end
        figsDataTrends{1+indROI} = subplot(3,2,indROI);
    else
        figsDataTrends{indROI} = figure(figBaseNum+indROI);
    end
    %% get data and trendline
    data = dataROIs(:,indROI) - desMatRem*betasRem(:,indROI);
    trendline = trendX*betaOI(indROI);
    
    %% the actual plotting
    plot(regColOI,data,'kx','MarkerSize',12); hold('on');
    plot(trendX,trendline,'r-','LineWidth',4);
    if(asOneFig)
        if(indROI==5||indROI==6)
            xlabel(colNameOI)
        else
            xticks([]);
        end
    else
        xlabel(colNameOI)
        title({[roiLabels{indROI},': "',colNameOI,'"-Trend over Residuals'];['(dataROI - \{',remStr,'\})']});
    end
    ylabel(roiLabels{indROI})
    %axis square
    %axis equal
    xlim([min([-2.875; min(trendX)*1.1]); max([1.375; max(trendX)*1.1])]); %xlim([-2.875; 1.375]); %
    %figsDataTrends{1+indROI}.PlotBoxAspectRatio = [5 2 1]; %[11 4 1]; %[16 9 1]; %[4 2 1];
end

end

%% plotFitResults
function figsRegResults = plotFitResults(fitResults,figNums,titleStr,plotMode,plotOnly,plotLevel)
% plot fitting results
nCols = length(fitResults.colNames_desMat);
if(~exist('figNums','var')||isempty(figNums))
    figNums = 420+(1:(nCols+2));
end
if(length(figNums)<(nCols+2))
    if(length(figNums)==1)
        figNums = figNums+(1:(nCols+2));
    else
        error(['figNums must be one or ',num2str(nCols+2),'==(nCols+2) numbers (integers).']);
    end
elseif(length(figNums)>(nCols+2))
    error(['figNums must be one or ',num2str(nCols+2),'==(nCols+2) numbers (integers).']);
end
    
if(~exist('titleStr','var')||isempty(titleStr))
    titleStr = '';
end

if(~exist('plotMode','var')||isempty(plotMode))
    plotMode = 'CI'; %plot data with CI (other option is SE ie standard errors
end

if(~exist('plotOnly','var')||isempty(plotOnly))
    plotOnly = 'ALL'; %plot all columns
end

%% get data for plotting
roiNames = fitResults.roiNames;
roiLabels = regexprep(roiNames,'_','\\_');
nROIs = length(roiNames);

allRsq = fitResults.stats(1,:);

betas = fitResults.betas;
betas_yLabels = fitResults.colNames_desMat(:);

pVals = [fitResults.pVals;fitResults.stats(3,:)];
pVal_yLabels = fitResults.colNames_desMat(:);
pVal_yLabels{end+1} = 'total';


%% do the plotting
figsRegResults = cell(3,1);
%% Rsquared
if(strcmp(plotLevel,'all'))
figsRegResults{1} = figure(figNums(1)); clf;
ax1 = subplot(1,1,1); 
bar(1:nROIs,allRsq.*100); hold on
xlabel('ROIs'); xticks(1:nROIs); xticklabels(roiLabels); xtickangle(-7.5);
ylabel('% explained variance');
ax1.FontSize = 24;
ylim([0, max(allRsq.*100)+1])
title(['R^2 of model "',titleStr,'"']);
end

%% betas and CI or SE
if(isfield(fitResults,'slopesComparison'))
    slopePairsOI = fitResults.slopesComparison.slopePairsOI;
    pValsPairs   = fitResults.slopesComparison.pVals;
else
    slopePairsOI = [];
    pValsPairs   = [];
end
switch(plotMode)
    case 'CI'
        %% betas with betaIs for p"Angle"
        for indCol = 1:size(fitResults.desMat,2) %length(fitResults.colNames_desMat)
            if(strcmp(plotOnly,'ALL')||strcmp(plotOnly,fitResults.colNames_desMat{indCol}))
                figsRegResults{1+indCol} = figure(figNums(1+indCol)); clf;
                plotBetasCIs(roiLabels,fitResults.colNames_desMat{indCol},fitResults.betas(indCol,:),squeeze(fitResults.betaIs(:,indCol,:)),fitResults.pVals(indCol,:),slopePairsOI,pValsPairs(indCol,:));
            end
        end
    case 'SE'
        %% betas with SEs for p"Angle"
        for indCol = 1:size(fitResults.desMat,2) %length(fitResults.colNames_desMat)
            if(strcmp(plotOnly,'ALL')||strcmp(plotOnly,fitResults.colNames_desMat{indCol}))
                figsRegResults{1+indCol} = figure(figNums(1+indCol)); clf;
                plotBetasSEs(roiLabels,fitResults.colNames_desMat{indCol},fitResults.betas(indCol,:),fitResults.stdErr(indCol,:),fitResults.pVals(indCol,:),slopePairsOI,pValsPairs(indCol,:));
            end
        end
    case '2SE'
        %% betas with 2*SE like betaIs for p"Angle"
        for indCol = 1:size(fitResults.desMat,2) %length(fitResults.colNames_desMat)
            if(strcmp(plotOnly,'ALL')||strcmp(plotOnly,fitResults.colNames_desMat{indCol}))
                figsRegResults{1+indCol} = figure(figNums(1+indCol)); clf;
                betas2SE = [fitResults.betas(indCol,:)-fitResults.stdErr(indCol,:); fitResults.betas(indCol,:)+fitResults.stdErr(indCol,:)];
                plotBetasCIs(roiLabels,fitResults.colNames_desMat{indCol},fitResults.betas(indCol,:),betas2SE,fitResults.pVals(indCol,:),slopePairsOI,pValsPairs(indCol,:),'StdErr');
            end
        end
    otherwise
        error(['ERROR: Unknown plot mode "',plotMode,'".']);
end

%% betas & p-Values
if(strcmp(plotLevel,'all')||strcmp(plotLevel,'sparse'))
figsRegResults{end} = figure(figNums(end)); clf;
subplot(1,2,1);
imagesc(betas); title('betas');
xticks(1:nROIs); xticklabels(roiLabels); xtickangle(-7.5);
yticks(1:length(betas_yLabels)); yticklabels(betas_yLabels);

subplot(1,2,2);
imagesc(pVals, [0 0.125]); title('p-Values');
xticks(1:nROIs); xticklabels(roiLabels); xtickangle(-7.5);
yticks(1:length(pVal_yLabels)); yticklabels(pVal_yLabels);

if(~isempty(titleStr))
    sgtitle(titleStr);
end
end

end

%% plotBetasSEs
function [ax] = plotBetasSEs(roiLabels,colNameOI,allBetas,allSEs,allPvals,slopePairsOI,pValsPairs)
% plot betas with SE and info about significance and differences of slopes.
if((~exist('slopePairsOI','var'))||(~exist('pValsPairs','var'))||isempty(slopePairsOI)||isempty(pValsPairs))
    slopePairsOI = [];
    pValsPairs   = [];
end

nROIs = length(roiLabels);

yLimits = [min(allBetas(:)-allSEs(:)); max(allBetas(:)+allSEs(:))];
ax = subplot(1,1,1); 
bar(1:nROIs,allBetas); hold on
for indROI = 1:nROIs
    if(allBetas(indROI)>=0)
        errorbar(indROI,allBetas(indROI),0,allSEs(indROI),'.','Color','r','LineWidth',4,'CapSize',24);
    else
        errorbar(indROI,allBetas(indROI),allSEs(indROI),0,'.','Color','r','LineWidth',4,'CapSize',24);
    end
end
for indROI = 1:nROIs
    %symbols formating
    currPval = allPvals(indROI);
    if(currPval<=0.001)
        symbol = '***';
        fontWeight = 'bold';
        fontSize = 16;
    elseif(currPval<=(0.05/nROIs)) %old <=0.005 0.05/nROIs==0.0083333...
        symbol = '**';
        fontWeight = 'bold';
        fontSize = 16;
    elseif(currPval<=0.01)
        symbol = '*';
        fontWeight = 'bold';
        fontSize = 14;
    elseif(currPval<=0.05)
        symbol = '+';
        fontWeight = 'normal';
        fontSize = 12;
    else
        symbol = '';
    end
    
    %write text and format pValue strings
    if(isempty(symbol))
        continue;
    else
        if(currPval<=0.001)
            pValStr = num2str(currPval,1);
        elseif(currPval<=0.01)
            pValStr = num2str(currPval,2);
        else
            pValStr = num2str(currPval,3);
        end
        
        if(allBetas(indROI)>=0)
            yPosSymbol  = allBetas(indROI)+allSEs(indROI);
            yPosPvalStr = yPosSymbol.*0.97; %allBetasI(2,indROI).*0.9;
            
            %update limits
            yLimits(2) = max([yLimits(2); yPosSymbol; yPosPvalStr]);
        else
            yPosSymbol  = allBetas(indROI)-allSEs(indROI);
            yPosPvalStr = yPosSymbol.*1.05; %allBetasI(1,indROI).*1.1;
            
            %update limits
            yLimits(1) = min([yLimits(1); yPosSymbol; yPosPvalStr]);
        end
        text(indROI+0.1, yPosSymbol,  symbol, 'FontWeight',fontWeight, 'FontSize',fontSize);
        text(indROI+0.1, yPosPvalStr,pValStr, 'FontWeight',  'normal', 'FontSize',12);
    end
end
ax.FontSize = 24;
xlabel('ROIs'); xticks(1:nROIs); xticklabels(roiLabels); xtickangle(-7.5);
ylabel('betas');
title(['betas with +SE of "',colNameOI,'"']);

if(~isempty(slopePairsOI)) %(isfield(fitResults,'slopesComparison'))
    for indTest = 1:size(slopePairsOI,1)
        currPair = slopePairsOI(indTest,:);
        currPVal = pValsPairs(indTest);
        if(currPVal<0.05)
            if(currPVal<=0.001)
                pValStr = num2str(currPVal,1);
            elseif(currPVal<=0.01)
                pValStr = num2str(currPVal,2);
            else
                pValStr = num2str(currPVal,3);
            end
            if((allBetas(currPair(1))>0)&&(allBetas(currPair(2))>0))
                yPosSymbol  = max([allSEs(currPair(1)); allSEs(currPair(2))]).*1.1;
                
                %update limits
                yLimits(2) = max([yLimits(2); yPosSymbol]);
            elseif((allBetas(currPair(1))<0)&&(allBetas(currPair(2))<0))
                yPosSymbol  = max([allSEs(currPair(1)); allSEs(currPair(2))]).*1.1;
                
                %update limits
                yLimits(1) = min([yLimits(1); yPosSymbol]);
            else
                yPosSymbol  = max(abs([allSEs(currPair(1)); allSEs(currPair(2))])).*1.1;
            end
            text((currPair(1)+currPair(2))/2, yPosSymbol,  ['|----(',pValStr,')----|'], 'FontWeight','bold', 'FontSize',16, 'HorizontalAlignment','center');
        end
    end
end
% set limits
if(yLimits(1)==yLimits(2))
    if(all(yLimits==0))
        yLimits = [0; 1];
    else
        yLimits = [yLimits(1) 2*yLimits(1)];
    end
end
yLimits(1) = min([-0.8; yLimits(1)]);
yLimits(2) = max([+0.9; yLimits(2)]);
ylim(yLimits); %ylim(yLimits.*1.1)
end

%% plotBetasCIs
function [ax] = plotBetasCIs(roiLabels,colNameOI,allBetas,allBetasI,allPvals,slopePairsOI,pValsPairs,intervalStr)
% plot betas with 95%-CI or 2*SE and info about significance and differences of slopes.
if((~exist('slopePairsOI','var'))||(~exist('pValsPairs','var'))||isempty(slopePairsOI)||isempty(pValsPairs))
    slopePairsOI = [];
    pValsPairs   = [];
end
if(~exist('intervalStr','var')||isempty(intervalStr))
    intervalStr = '95%-CI';
end

nROIs = length(roiLabels);

yLimits = [min(allBetasI(:)); max(allBetasI(:))];
ax = subplot(1,1,1); 
bar(1:nROIs,allBetas); hold on
errorbar(1:nROIs,allBetas,allBetas-allBetasI(1,:),allBetas-allBetasI(2,:),'.','Color','r','LineWidth',4,'CapSize',24);
for indROI = 1:nROIs
    %symbols formating
    currPval = allPvals(indROI);
    if(currPval<=0.001)
        symbol = '***';
        fontWeight = 'bold';
        fontSize = 16;
    elseif(currPval<=(0.05/nROIs)) %old <=0.005 0.05/nROIs==0.0083333...
        symbol = '**';
        fontWeight = 'bold';
        fontSize = 16;
    elseif(currPval<=0.01)
        symbol = '*';
        fontWeight = 'bold';
        fontSize = 14;
    elseif(currPval<=0.05)
        symbol = '+';
        fontWeight = 'normal';
        fontSize = 12;
    else
        symbol = '';
    end
    
    %write text and format pValue strings
    if(isempty(symbol))
        continue;
    else
        if(currPval<=0.001)
            pValStr = num2str(currPval,1);
        elseif(currPval<=0.01)
            pValStr = num2str(currPval,2);
        else
            pValStr = num2str(currPval,3);
        end
        
        if(allBetas(indROI)>0)
            yPosSymbol  = allBetasI(2,indROI);
            yPosPvalStr = yPosSymbol.*0.97; %allBetasI(2,indROI).*0.9;
            
            %update limits
            yLimits(2) = max([yLimits(2); yPosSymbol; yPosPvalStr]);
        else
            yPosSymbol  = allBetasI(1,indROI);
            yPosPvalStr = yPosSymbol.*1.05; %allBetasI(1,indROI).*1.1;
            
            %update limits
            yLimits(1) = min([yLimits(1); yPosSymbol; yPosPvalStr]);
        end
        text(indROI+0.1, yPosSymbol,  symbol, 'FontWeight',fontWeight, 'FontSize',fontSize);
        text(indROI+0.1, yPosPvalStr,pValStr, 'FontWeight',  'normal', 'FontSize',12);
    end
end
ax.FontSize = 24;
xlabel('ROIs'); xticks(1:nROIs); xticklabels(roiLabels); xtickangle(-7.5);
ylabel('betas');
title(['betas with ',intervalStr,' of "',colNameOI,'"']);

if(~isempty(slopePairsOI)) %(isfield(fitResults,'slopesComparison'))
    for indTest = 1:size(slopePairsOI,1)
        currPair = slopePairsOI(indTest,:);
        currPVal = pValsPairs(indTest);
        if(currPVal<0.05)
            if(currPVal<=0.001)
                pValStr = num2str(currPVal,1);
            elseif(currPVal<=0.01)
                pValStr = num2str(currPVal,2);
            else
                pValStr = num2str(currPVal,3);
            end
            if((allBetas(currPair(1))>0)&&(allBetas(currPair(2))>0))
                yPosSymbol  = max([allBetasI(2,currPair(1)); allBetasI(2,currPair(2))]).*1.1;
                
                %update limits
                yLimits(2) = max([yLimits(2); yPosSymbol]);
            elseif((allBetas(currPair(1))<0)&&(allBetas(currPair(2))<0))
                yPosSymbol  = max([allBetasI(1,currPair(1)); allBetasI(1,currPair(2))]).*1.1;
                
                %update limits
                yLimits(1) = min([yLimits(1); yPosSymbol]);
            else
                yPosSymbol  = max(abs([allBetasI(2,currPair(1)); allBetasI(2,currPair(2))])).*1.1;
            end
            text((currPair(1)+currPair(2))/2, yPosSymbol,  ['|----(',pValStr,')----|'], 'FontWeight','bold', 'FontSize',16, 'HorizontalAlignment','center');
        end
    end
end
% set limits
if(yLimits(1)==yLimits(2))
    if(all(yLimits==0))
        yLimits = [0; 1];
    else
        yLimits = [yLimits(1) 2*yLimits(1)];
    end
end
yLimits(1) = min([-0.8; yLimits(1)]);
yLimits(2) = max([+1;   yLimits(2)]);
ylim(yLimits); %ylim([-0.8; 1]); %ylim(yLimits.*1.1);
end

%% fitAllROIs
function fitResults = fitAllROIs(dataROIs,roiNames,regressorsTable,desMatSettings)
% fit model to all ROIs
if(~exist('desMatSettings','var')||isempty(desMatSettings))
    desMatSettingsOrg = {'pMVS','Sex*Handedness','Age','dt_RSfMRI','const'};
else
    desMatSettingsOrg = desMatSettings;
end

%% make sure desMatSettings does not contain _ALL suffix
desMatSettings = replaceALLtypes(desMatSettingsOrg);

%% generate design matrix
[desMat, colNames_desMat, idxAngleCol] = genDesMat(regressorsTable,desMatSettings);
nCols = size(desMat,2); %length(colNames_desMat);

%% check if desMat is rank deficient (only output warning once not at every ROI regression, as the design matrix is anyways the same)
% [hasFullRank,xrank] = checkDesMatHasFullRank(desMat);
% if(~hasFullRank)
%     if(ismember('useFixRankDefReg',evalin('base','who')))
%         useFixRankDefReg = evalin('base','useFixRankDefReg');
%     else
%         useFixRankDefReg = true;
%     end
% %     if(useFixRankDefReg)
% %         warning('fitAllROIs:desMatRankDeficient',...
% %        ['WARNING: it appears that the design matrix is rank deficient nCols==',num2str(nCols),' but rank(desMat)==',num2str(xrank),'! Will try to compensate for rank deficiency.']);
% %     else
% %         warning('fitAllROIs:desMatRankDeficient',...
% %        ['WARNING: it appears that the design matrix is rank deficient nCols==',num2str(nCols),' but rank(desMat)==',num2str(xrank),'!']);
% %     end
% end

%% fit each ROI
nROIs = size(dataROIs,2);
betas = zeros(  nCols,nROIs);
stdErr= zeros(  nCols,nROIs);
rSq   = zeros(        nROIs,1);
AIC   = zeros(        nROIs,1);
BIC   = zeros(        nROIs,1);
betaIs= zeros(2,nCols,nROIs);
pVals = ones(   nCols,nROIs);
stats = nan(        6,nROIs);
stats_rob = cell(     nROIs,1);
for indROI = 1:nROIs
    [betas(:,indROI),betaIs(:,:,indROI),pVals(:,indROI),stats(:,indROI),stats_rob{indROI}] = altRobustFit(dataROIs(:,indROI),desMat);
    rSq(     indROI) = stats(1,indROI);
    AIC(     indROI) = stats(5,indROI); 
    BIC(     indROI) = stats(6,indROI);
    stdErr(:,indROI) = stats_rob{indROI}.se;
end

%% assign output struct
fitResults = struct('settings',struct('regressorsTable',regressorsTable,...
                                      'desMatSettings',{desMatSettings},...
                                      'desMatSettingsOrg',{desMatSettingsOrg}),...
                    'desMat',desMat,...
                    'colNames_desMat',{colNames_desMat},...
                    'idxAngleCol',idxAngleCol,...
                    'dataROIs',dataROIs,...
                    'roiNames',{roiNames},...
                    'betas',betas,...
                    'pVals',pVals,...
                    'stdErr',stdErr,...
                    'AIC',AIC,...
                    'BIC',BIC,...
                    'rSq',rSq,...
                    'stats',stats,...
                    'betaIs',betaIs,...
                    'stats_rob',{stats_rob});
 
%% compare slopes
fitResults = compareSlopesLR(fitResults); % do slope comparison --> plot between L & R areas should be done by plotFitResults

end

%% compareSlopesLR
function fitResults = compareSlopesLR(fitResults)
% compare the slopes between the ROIs for left and right sides
nSubjs = size(fitResults.dataROIs,1);

slopePairsOI = getSlopesOI(fitResults.roiNames);
%colNames_desMat = fitResults.colNames_desMat;
pVals = nan(  size(fitResults.desMat,2),size(slopePairsOI,1)); %nan(  length(colNames_desMat),size(slopePairsOI,1));
tVals = zeros(size(fitResults.desMat,2),size(slopePairsOI,1)); %zeros(length(colNames_desMat),size(slopePairsOI,1));
dfe   = zeros(size(fitResults.desMat,2),size(slopePairsOI,1)); %zeros(length(colNames_desMat),size(slopePairsOI,1));
for indTest = 1:size(slopePairsOI,1)
    slopeIdxL = slopePairsOI(indTest,1);
    slopeIdxR = slopePairsOI(indTest,2);
    for indCol = 1:size(fitResults.desMat,2) %length(colNames_desMat)
        beta_1 = fitResults.betas(indCol,slopeIdxL);
        beta_2 = fitResults.betas(indCol,slopeIdxR);
        
        se_1   = fitResults.stats_rob{slopeIdxR}.se(indCol);
        se_2   = fitResults.stats_rob{slopeIdxR}.se(indCol);
        [pVals(indCol,indTest),tVals(indCol,indTest),dfe(indCol,indTest)] = compareTwoSlopes(beta_1,se_1,beta_2,se_2,nSubjs);
    end
end

fitResults.slopesComparison = struct('slopePairsOI',slopePairsOI,...
                                     'pVals',pVals,...
                                     'tVals',tVals,...
                                     'dfe',dfe);
end

%% getSlopesOI
function slopePairsOI = getSlopesOI(roiNames)
% find names which end in left or right (or L & R)
% remove/ignore those that do not have those endings
% for the remaining, remove the left/right string and compare the strings
% --> this should yield pairs
% For each pair determine which is left --> this is the first entry the other
% one must be right --> 2nd entry
% return the nTests-x-2 matrix of indices for testing.

idxLeft  = find(cellfun(@(str) ~isempty(str),regexp(roiNames,' left| Left')));
idxRight = find(cellfun(@(str) ~isempty(str),regexp(roiNames,' right| Right')));


slopePairsOI = matchIdxLR(regexprep(roiNames,' left| Left| right| Right',''),idxLeft,idxRight);
end

%% matchIdxLR
function slopePairsOI = matchIdxLR(roiNamesOI,idxLeft,idxRight)
% compare the names of the ROIs for left and right
slopePairsOI = zeros(length(idxLeft),2);
for indPair = 1:length(idxLeft)
    if(strcmp(roiNamesOI{idxLeft(indPair)},roiNamesOI{idxRight(indPair)}))
        slopePairsOI(indPair,:) = [idxLeft(indPair), idxRight(indPair)];
    else
        keyboard;
    end
end
end

%% compareTwoSlopes
function [pVal,tVal,dfe] = compareTwoSlopes(beta_1,se_1,beta_2,se_2,nTotal)
% do comparison of slopes
% see appendix section 2.2. of
% Cohen, J., Cohen, P., West, S.G., and Aiken, L.S. (2003). 
% Applied Multiple Regression/Correlation Analysis for the Behavioral Sciences (3rd edition). 
% Mahwah, NJ: Lawrence Earlbaum Associates.

dB = beta_2 - beta_1; 
sSq= sqrt(se_1^2 + se_2^2); 
tVal = dB/sSq;
dfe = nTotal - 4; %two slopes and two std-errors are estimated --> 4 df are gone
pVal = 2 * tcdf(-abs(tVal), dfe);

end

%% replaceALLtypes
function newDesMatSettings = replaceALLtypes(desMatSettings)
% in case there is 'Age_ALL' or 'dt_RSfMRI_ALL', then replace this with
% zscore(Age) AND zscore(tiedrank(Age)) and analog for 'dt_RSfMRI'

selAgeAll = cellfun(@(str) strcmp(str,'Age_ALL'),desMatSettings);
hasAgeAll = any(selAgeAll);
sel_dt_RSfMRI_ALL = cellfun(@(str) strcmp(str,'dt_RSfMRI_ALL'),desMatSettings);
has_dt_RSfMRI_ALL = any(sel_dt_RSfMRI_ALL);
if(hasAgeAll&&has_dt_RSfMRI_ALL)
    selAge_dt_RSfMRI = selAgeAll|sel_dt_RSfMRI_ALL;
    selOthers = ~selAge_dt_RSfMRI;
    if(any(selOthers))
        newDesMatSettings = desMatSettings(selOthers); %keep all other regressors
        newDesMatSettings{end+1} = 'tiedrank(Age)';
        newDesMatSettings{end+1} = 'Age';
        newDesMatSettings{end+1} = 'tiedrank(dt_RSfMRI)';
        newDesMatSettings{end+1} = 'dt_RSfMRI';
    else %there are no other regressors besides Age & dt_RSfMRI
        newDesMatSettings = {'tiedrank(Age)','demean(Age)','tiedrank(dt_RSfMRI)','demean(dt_RSfMRI)'};
    end
elseif(~hasAgeAll&&has_dt_RSfMRI_ALL) %only dt_RSfMRI_ALL
    selOthers = ~sel_dt_RSfMRI_ALL;
    if(any(selOthers))
        newDesMatSettings = desMatSettings(selOthers); %keep all other regressors
        newDesMatSettings{end+1} = 'tiedrank(dt_RSfMRI)';
        newDesMatSettings{end+1} = 'dt_RSfMRI';
    else %there are no other regressors besides dt_RSfMRI
        newDesMatSettings = {'tiedrank(dt_RSfMRI)','demean(dt_RSfMRI)'};
    end
elseif(hasAgeAll&&~has_dt_RSfMRI_ALL) %only Age_ALL
    selOthers = ~selAgeAll;
    if(any(selOthers))
        newDesMatSettings = desMatSettings(selOthers); %keep all other regressors
        newDesMatSettings{end+1} = 'tiedrank(Age)';
        newDesMatSettings{end+1} = 'Age';
    else %there are no other regressors besides Age 
        newDesMatSettings = {'tiedrank(Age)','Age'};
    end
else
    newDesMatSettings = desMatSettings; %leave unchanged
end
end

%% genDesMat
function [desMat, colNames, idxAngleCol] = genDesMat(regressorsTable,regNames)
% generate design matrix, if a regressor is continuous type, then zscore in
% order to get normalized betas in the regression later.
% "Constant"/Categorical regressors are output as multiple columns of zeros and ones.
% 
%outputs:
%        desMat      <-- the design matrix for robust regression
%        colNames    <-- column names (for plotting and as a check for the correctness of the assignment)
%        idxAngleCol <-- the column index for the "angle" regressor (pMVS or pReids).
%                        There should be only one but mainly this is just so we
%                        can use this to extract the correct column for plotting
%                        the data and trendlines.

%% inputs
assert(istable(regressorsTable),'ERROR: regressorsTable must be table!');
assert(iscellstr(regNames),'ERROR: regNames must be a cellstring!');

%% fill design matrix
desMat = []; %init and fill up
colNames = {}; %init and fill up
idxAngleCol = []; %init and fill up
for indR = 1:length(regNames)
    [currCols, currColNames, isAngleReg] = getColumnsForDesMat(regressorsTable,regNames{indR});
    if(isAngleReg)
        idxAngleCol = [idxAngleCol; size(desMat,2)+(1:size(currCols,2))'];
    end
    desMat  = [desMat, currCols];
    colNames= [colNames; currColNames(:)];
end

%% orthogonalize
if(~isempty(idxAngleCol))
    desMat = orthogonalize(desMat,idxAngleCol);
end

end

%% getColumnsForDesMat
function [colsData, colNames, isAngleReg] = getColumnsForDesMat(regressorsTable,regName)
% pick regressor, 
% if continuous use zscore
%    check if it is one of the Angle cols --> isAngleCol==true
% else make 0 & 1 regressors.
%
%

%% inputs
assert(ischar(regName),'ERROR: regName must be a char-vector!');

%% find appropriate output for regName
switch(regName)
    case {'pReidsL';'pReidsR';'pReidsAv';'pMVSL';'pMVSR';'pMVS'}
        currData = regressorsTable.(regName);
        combOutliers = sum([isoutlier(currData,'mean'),isoutlier(currData,'median'),isoutlier(currData,'quartiles'),isoutlier(currData,'grubbs'),isoutlier(currData,'gesd')],2);
        if(any(combOutliers>2))
            currData(combOutliers>2) = nan;
        end
        colsData = nanzscore(currData);
        colNames = {regName};
        isAngleReg = true;
    case 'Sex'
        colsData = [(regressorsTable.Sex=='female'),...
                    (regressorsTable.Sex=='male')];
        colNames = {'female';'male'};
        isAngleReg = false;
    case 'SexDiff'
        colsData = (regressorsTable.Sex=='female')-(regressorsTable.Sex=='male');
        colNames = {'[female-male]'};
        isAngleReg = false;
    case 'Levels(Sex)'
        colsData = (regressorsTable.Sex=='female')+2.*(regressorsTable.Sex=='male');
        colNames = {'[1==Female;2==Male]'};
        isAngleReg = false;
    case 'Sex*Handedness'
        colsData = [(regressorsTable.Sex=='female').*(regressorsTable.Handedness=='left'), ...
                    (regressorsTable.Sex==  'male').*(regressorsTable.Handedness=='left'), ...
                    (regressorsTable.Sex=='female').*(regressorsTable.Handedness=='right'), ...
                    (regressorsTable.Sex==  'male').*(regressorsTable.Handedness=='right')];
        colNames = {'Female-LH';'Male-LH';'Female-RH';'Male-RH'};
        isAngleReg = false;
    case 'Levels(Sex*Handedness)'
        colsData = [(regressorsTable.Sex=='female').*(regressorsTable.Handedness=='left'), ...
                    (regressorsTable.Sex==  'male').*(regressorsTable.Handedness=='left'), ...
                    (regressorsTable.Sex=='female').*(regressorsTable.Handedness=='right'), ...
                    (regressorsTable.Sex==  'male').*(regressorsTable.Handedness=='right')];
        colsData = sum(colsData.*[1,2,3,4],2);
        colNames = {'[1==Female-LH;2==Male-LH;3==Female-RH;3==Male-RH]'};
        isAngleReg = false;
    case 'Handedness'
        colsData = [(regressorsTable.Handedness=='right'),...
                    (regressorsTable.Handedness=='left')];
        colNames = {'RH';'LH'};
        isAngleReg = false;
    case 'HandednessDiff'
        colsData = (regressorsTable.Handedness=='right')-(regressorsTable.Handedness=='left');
        colNames = {'RH-LH'};
        isAngleReg = false;
    case 'Levels(Handedness)'
        colsData = (regressorsTable.Handedness=='right') + 2.*(regressorsTable.Handedness=='left');
        colNames = {'[1==RH;2==LH]'};
        isAngleReg = false;
    case {'Age','demean(Age)','zscore(Age)'}
        colsData = zscore(regressorsTable.Age);
        colNames = {'Age'};
        isAngleReg = false;
    case 'tiedrank(Age)'
        colsData = zscore(tiedrank(regressorsTable.Age));
        colNames = {'tiedrank(Age)'};
        isAngleReg = false;
    case {'dt_RSfMRI','demean(dt_RSfMRI)','zscore(dt_RSfMRI)'}
        colsData = zscore(regressorsTable.dt_RSfMRI);
        colNames = {'dt-RSfMRI'};
        isAngleReg = false;
    case 'tiedrank(dt_RSfMRI)'
        colsData = zscore(tiedrank(regressorsTable.dt_RSfMRI));
        colNames = {'tiedrank(dt-RSfMRI)'};
        isAngleReg = false;
    case 'medianSplit(dt_RSfMRI)'
        dtMedian = median(regressorsTable.dt_RSfMRI);
        colsData = [double(regressorsTable.dt_RSfMRI<dtMedian) + 2.*double(regressorsTable.dt_RSfMRI>=dtMedian)];
        colNames = {['dt-RSfMRI<',num2str(dtMedian),';dt-RSfMRI>=',num2str(dtMedian)]};
        isAngleReg = false;
    case 'group2(dt_RSfMRI)'
        colsData = [double(regressorsTable.dt_RSfMRI<34) + 2.*double(regressorsTable.dt_RSfMRI>=34)];
        colNames = {'dt-RSfMRI<34;dt-RSfMRI>=34'};
        isAngleReg = false;
    case 'group3(dt_RSfMRI)'
        colsData = [double(regressorsTable.dt_RSfMRI<10) + 2.*double(regressorsTable.dt_RSfMRI>=10&regressorsTable.dt_RSfMRI<=20) + 3.*double(regressorsTable.dt_RSfMRI>20)];
        colNames = {'dt-RSfMRI<10;dt-RSfMRI>=10&<=20;dt-RSfMRI>20'};
        isAngleReg = false;
    case 'const'
        colsData = ones(size(regressorsTable,1),1);
        colNames = {'const'};
        isAngleReg = false;
    otherwise
        error(['Unknown regressor type "',regName,'"! Please check your inputs.']);
end
end

%% orthogonalize
function newDesMat = orthogonalize(desMat,idxAngleCol)
% orthogonalize without damaging the angle columns or categorical/constant(ish)/0&1 columns
% use robust regression to remove contribution of angle columns

nCols = size(desMat,2);
newDesMat = zeros(size(desMat)); %init
pAngle = desMat(:,idxAngleCol);
idxCat = findCatCols(desMat);
for ind = 1:nCols   
    if(any(idxAngleCol==ind)||any(idxCat==ind)) %don't damage
        newDesMat(:,ind) = desMat(:,ind);
    else
        [~,~,~,~,~,res] = altRobustFit(desMat(:,ind),pAngle);
        newDesMat(:,ind) = nanzscore(res);
    end
end

end

%% findCatCols
function idxCat = findCatCols(desMat)
% find the columns that have only zeros and ones in them
% if none are present return ZERO

idxCat = []; %init and fill up
for ind = 1:size(desMat,2)
    uniqueVals = unique(desMat(:,ind));
    uniqueVals(uniqueVals==0|uniqueVals==1) = [];
    if(isempty(uniqueVals))
        idxCat = [idxCat; ind];
    end
end
if(isempty(idxCat))
    idxCat = 0;
end
end

%% altRobustFit
function [betas,betaIs,pVals,stats,stats_rob,res,resI] = altRobustFit(y,X,wfun,const)
% This is an alternative version of the robustfit.m function.
% It uses robustfit with const=='off' as it is assumed that the design matrix X,
% already has columns with "1"s.
%
% The outputs are adjusted such that use is easier and additional infos like for 
% regress.m are output.
%
%NB(V2.0+): there is only a warning for rank deficiency displayed and results 
%           are adjusted if possible.  
%
%Usage:
%      [betas,betaIs,pVals,stats,stats_rob,res,resI] = altRobustFit(y,X,wfun,const)
%      [betas,betaIs,pVals,stats,stats_rob] = altRobustFit(y,X,'bisquare','off');
%      [betas,betaIs,pVals,stats,stats_rob] = altRobustFit(y,X); %same as above
%
%Author: Rainer.Boegle@googlemail.com
%V3.0: 20.05.2020
%Comment: V3.0(20.05.2020): added AIC & BIC estimation. V2.0(22.04.2020): added behavior for call with hasFullRank. V1.0(10.04.2020): initial implementation.

%% settings
if(~exist('const','var')||isempty(const))
    const = 'off';
end
if(~strcmp(const,'off'))
    error('ERROR: this function cannot work if constant is added! (Sorry)');
end

if(~exist('wfun','var')||isempty(wfun))
    wfun = 'huber'; %'bisquare'; %'andrews'; %'fair'; %'cauchy'; %'welsch'; %'talwar'; %'logistic'; %'ols'; %
end
switch(wfun)
    case 'bisquare' %(default)	w = (abs(r)<1) .* (1 - r.^2).^2	
        tune = 4.685;
    case 'andrews'	%w = (abs(r)<pi) .* sin(r) ./ r	
        tune = 1.339;
    case 'cauchy'	%w = 1 ./ (1 + r.^2)	
        tune = 2.385;
    case 'fair'     %w = 1 ./ (1 + abs(r))	
        tune = 1.400;
    case 'huber'	%w = 1 ./ max(1, abs(r))	
        tune = 1.345;
    case 'logistic'	%w = tanh(r) ./ r	
        tune = 1.205;
    case 'ols'      %Ordinary least squares (no weighting function)	None
        tune = [];
    case 'talwar'	%w = 1 * (abs(r)<1)	
        tune = 2.795;
    case 'welsch'	%w = exp(-(r.^2))	
        tune = 2.985;
end

%% check rank of design matrix
[hasFullRank,xrank] = checkDesMatHasFullRank(X); %hasFullRank = true; %for debug
if(ismember('useFixRankDefReg',evalin('base','who')))
    useFixRankDefReg = evalin('base','useFixRankDefReg');
else
    useFixRankDefReg = true;
end

%% run robust regression
priorW = ones(size(y));
doWarn = false;
[b_rob, stats_rob] = robustfit(X,y,wfun,tune,const,priorW,doWarn);
if(~hasFullRank&&useFixRankDefReg) %try to compensate for rank deficiency
    [b_rob, stats_rob] = fixRankDefRegression(b_rob,stats_rob,X,wfun,tune,const);
end

%% get R-square
sse     = stats_rob.dfe * stats_rob.robust_s^2; %sum of squares error
if(any(isnan(X(:))))
    phat = zeros(size(X,1),1);
    for indS = 1:size(X,1)
        phat(indS) = sum(X(indS,:)'.*b_rob,'omitnan');
    end
else
    phat    = X*b_rob; %prediction
end
ssr     = sum((phat-mean(phat,'omitnan')).^2,'omitnan'); %sum of squares regression %old: norm(phat-mean(phat))^2;
tss     = sse + ssr; %ovserved sum of squares
rsquare = 1 - sse/tss; %https://de.mathworks.com/matlabcentral/answers/93865-how-do-i-compute-the-r-square-statistic-for-robustfit-using-statistics-toolbox-7-0-r2008b

%% estimate AIC & BIC (I am more comfortable with using only the number of subjects that were non NaN and the rank of the design matrix after adjustment instead of the raw numbers)
nSubjsAll  = size(X,1);
nNaN       = max(sum(isnan(X),1),[],2);
nSubjs = nSubjsAll - nNaN;
nParams = xrank; %size(X,2);

%AIC = n ln(SSE) ? n ln(n) + 2p
AIC = nSubjs*log(sse) - nSubjs*log(nSubjs) + 2*nParams;

%BIC = n ln(SSE) - n ln(n) + ln(n)p
BIC = nSubjs*log(sse) - nSubjs*log(nSubjs) + log(nSubjs)*nParams;


%% assemble outputs
pVals = stats_rob.p(:);
betas = b_rob(:); %parameters
pCI  = 1-0.025; %95%-CI %1-0.1; %80&-CI %1-0.05; %90%-CI
tVal = tinv(pCI,stats_rob.dfe);
betaIs = [betas-(tVal(1).*stats_rob.se(:)), betas+(tVal(1).*stats_rob.se(:))]';
%betasI = [b_rob(:)-stats_rob.se(:), b_rob(:)+stats_rob.se(:)]; %confidence interval for parameters NB: there is a tvalue factor missing here, will fix this later.
if(any(isnan(X(:))))
    res = y - phat;
else
    res = stats_rob.resid(:); %residuals
end
nu = stats_rob.dfe; %max(0,length(y)-size(X,2));                % Residual degrees of freedom
if nu ~= 0
    tval = tinv((1-0.025),nu);
else
    tval = 0;
end 
ser    = (res(:)-mean(res(:),'omitnan')).^2;
resI   = [(res-tval.*ser) (res+tval.*ser)]; % Create confidence intervals for residuals.
s2     = stats_rob.robust_s^2;
if size(X,2) > 1
    F  = (ssr/(size(X,2)-1))/s2;      % F statistic for regression
else
    F  = NaN;
end
prob = fpval(F,size(X,2)-1,nu); % Significance probability for regression (the full model not the individual regressors)
stats  = [rsquare; F; prob; s2; AIC; BIC];

end

%% fpval
function p = fpval(x,df1,df2)
%FPVAL F distribution p-value function.
%   P = FPVAL(X,V1,V2) returns the upper tail of the F cumulative distribution
%   function with V1 and V2 degrees of freedom at the values in X.  If X is
%   the observed value of an F test statistic, then P is its p-value.
%
%   The size of P is the common size of the input arguments.  A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also FCDF, FINV.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.

%   Copyright 2010 The MathWorks, Inc. 


if nargin < 3, 
    error(message('stats:fpval:TooFewInputs')); 
end

xunder = 1./max(0,x);
xunder(isnan(x)) = NaN;
p = fcdf(xunder,df2,df1);

end

%% checkDesMatHasFullRank
function [hasFullRank,xrank] = checkDesMatHasFullRank(desMat)
% check if desMat is rank deficient
nSubj = size(desMat,1);
nCols = size(desMat,2);

if(any(isnan(desMat(:)))) %remove NaN if necessary
    [anybad wasnan X] = statremovenan(desMat);
    if (anybad==2)
        error(message('checkDesMatHasFullRank:InputSizeMismatch'));
    end
    % Find the least squares solution.
    [~,R] = qr(X,0);
else
    % Find the least squares solution.
    [~,R] = qr(desMat,0);
end

if isempty(R) % can only happen if p=0 as we know n>p>=0
    tol = 1;
else
    tol = abs(R(1)) * max(nSubj,nCols) * eps(class(R));
end
xrank = sum(abs(diag(R)) > tol);
if(xrank~=nCols)
    hasFullRank = false;
else
    hasFullRank = true;
end

end

%% statremovenan (private in MATLAB but I need it here)
function [badin,wasnan,varargout]=statremovenan(varargin)
%STATREMOVENAN Remove NaN values from inputs

%   Copyright 1993-2012 The MathWorks, Inc.


[badin,wasnan,varargout{1:nargout-2}] = internal.stats.removenan(varargin{:});
end

%% fixRankDefRegression
function [b_robNew, stats_robNew] = fixRankDefRegression(b_rob,stats_rob,X,wfun,tune,const)
% try to fix rank deficiency by identifying prolematic columns (usually only one)
% If only one column, then generate residuals and retry fitting problematic column on residuals.
% If MORE than one column, enter debug mode and give out warning.

res   = stats_rob.resid(:); %residuals
pVals = stats_rob.p(:); %p-values

idxColOI = find(isnan(pVals));
if(isempty(idxColOI))
    error('Somehow, the desMat is rank deficient but there are no pVals that are NaN??? This is strange, we should not have been in this execution path and not find any isNaN(pVals)??? Check what is going on here!');
elseif(length(idxColOI)==1)
    [b_ColOI, stats_ColOI] = robustfit(X(:,idxColOI),res,wfun,tune,const);

    b_robNew = b_rob;
    b_robNew(idxColOI) = b_ColOI;
    stats_robNew = stats_rob;
    stats_robNew.se(idxColOI) = stats_ColOI.se;
    stats_robNew.p(idxColOI)  = stats_ColOI.p;
    stats_robNew.t(idxColOI)  = stats_ColOI.t;
    %keyboard; %try to figure out what to do
elseif(length(idxColOI)==2)
    [b_ColOI, stats_ColOI] = robustfit(X(:,idxColOI),res,wfun,tune,const);

    b_robNew = b_rob;
    b_robNew(idxColOI) = b_ColOI;
    stats_robNew = stats_rob;
    stats_robNew.se(idxColOI) = stats_ColOI.se;
    stats_robNew.p(idxColOI)  = stats_ColOI.p;
    stats_robNew.t(idxColOI)  = stats_ColOI.t;
    %warning('fixRankDefRegression:NANcolsFound',...
    %        'WARNING: there are TWO columns with NaN pValues!!! Will treat it just as if there was 1, however it could be a problem!');
else
    disp('WARNING: there are multiple columns with NaN pValues!!! Check what is going on!');
    disp('WARNING: will not even try to fix results, -sorry.');
    b_robNew = b_rob;
    stats_robNew = stats_rob;
    keyboard; %default behavior should better be to just return and keep the problem...
    return;    
end

end

%% checkPlotCorrs
function [figsCorr,corrStruct] = checkPlotCorrs(roiNames,dataROIs,covarNames,covarsMat,figNum)
% correlate covars and ROI data, AND individually with each other.
if(~exist('figNum','var')||isempty(figNum))
    figNum = [1,2];
elseif((length(figNum)==1))
    figNum = [figNum, figNum+1];
elseif((length(figNum)>2))
    error('Figure number must be one input (expanded to two) or two.');
end

[rhoPearson_DataCovs,pValPearson_DataCovs] = corr(dataROIs,covarsMat);
[rhoPearson_Data    ,pValPearson_Data]     = corr(dataROIs);
[rhoPearson_Covs    ,pValPearson_Covs]     = corr(covarsMat);

[rhoSpearman_DataCovs,pValSpearman_DataCovs] = corr(dataROIs,covarsMat,'type','Spearman');
[rhoSpearman_Data    ,pValSpearman_Data]     = corr(dataROIs,'type','Spearman');
[rhoSpearman_Covs    ,pValSpearman_Covs]     = corr(covarsMat,'type','Spearman');


[rhoPearson_all, pValPearson_all]  = corr([dataROIs,covarsMat]); 
[rhoSpearman_all,pValSpearman_all] = corr([dataROIs,covarsMat],'type','Spearman');

%% assign struct
corrStruct = struct('Data_x_Cov',struct('Pearson', struct('rho',rhoPearson_DataCovs, 'pVals',pValPearson_DataCovs),...
                                        'Spearman',struct('rho',rhoSpearman_DataCovs,'pVals',pValSpearman_DataCovs)),...
                      'DataROIs',struct('Pearson', struct('rho',rhoPearson_Data,     'pVals',pValPearson_Data),...
                                        'Spearman',struct('rho',rhoSpearman_Data,    'pVals',pValSpearman_Data)),...
                        'Covars',struct('Pearson', struct('rho',rhoPearson_Covs,     'pVals',pValPearson_Covs),...
                                        'Spearman',struct('rho',rhoSpearman_Covs,    'pVals',pValSpearman_Covs)));
                                    
%% plot
allCols = [roiNames;covarNames];
figsCorr{2} = figure(figNum(1)); clf;
plotCorrs(rhoPearson_all,pValPearson_all,allCols,'Pearson')

figsCorr{2} = figure(figNum(2)); clf;
plotCorrs(rhoSpearman_all,pValSpearman_all,allCols,'Spearman')

end

%% plotCorrs
function plotCorrs(rho_all,pVal_all,allCols,nameStr)
imagesc(rho_all,[-1 1]); title([nameStr,'-Corr: dataROIs-x-Covars']); xticks(1:length(allCols));    xticklabels(allCols); yticks(1:length(allCols)); yticklabels(allCols);
for indC = 1:size(pVal_all,2)
    for indR = 1:size(pVal_all,1)
        currPval = pVal_all(indR,indC);
        if(currPval<=0.001)
            symbol = '***';
            fontSize = 16;
            fontWeight = 'bold';
        elseif(currPval<=0.005)
            symbol = '**';
            fontSize = 16;
            fontWeight = 'bold';
        elseif(currPval<=0.01)
            symbol = '*';
            fontSize = 14;
            fontWeight = 'bold';
        elseif(currPval<0.05)
            symbol = '+';
            fontSize = 12;
            fontWeight = 'normal';
        else
            symbol = '';
        end
        if(isempty(symbol))
            continue;
        else
            if(currPval<10^-3)
                pValStr = num2str(currPval,1);
            elseif(currPval<10^-2)
                pValStr = num2str(currPval,2);
            else
                pValStr = num2str(currPval,3);
            end
            pValYshift = 0.2;
            text(indR,indC,symbol,'FontWeight',fontWeight,'FontSize',fontSize,'HorizontalAlignment','center');
            text(indR,indC+pValYshift,pValStr,'FontWeight',fontWeight,'FontSize',fontSize,'HorizontalAlignment','center');
        end
    end
end
end

%% loadData
function dataReAnalysis = loadData(pathDataOI,drop2Oldest)
% load the data or load and drop the two oldest subjects


if(drop2Oldest)
    load(pathDataOI,'dataReAnalysis');
    Age = dataReAnalysis.regressorsTable.Age;
    uniqueAge = unique(Age);
    limitAge = uniqueAge(end-1);
    sel = (Age<limitAge);
    
    dataReAnalysis.regressorsTable = dataReAnalysis.regressorsTable(sel,:);
    dataReAnalysis.dataROIs        = dataReAnalysis.dataROIs(sel,:);
else
    load(pathDataOI,'dataReAnalysis');
end
end

%% prepend
function outCell = prepend(prependElements,inCell)
% prepend the prependElements

assert(iscell(inCell),'ERROR: inCell must be a cell-vector!');
assert(isvector(inCell),'ERROR: inCell must be a vector (cell-vector)!');

nIn  = length(inCell);
if(iscell(prependElements))
    assert(isvector(prependElements),'ERROR: prependElements must be a single element or a cell-vector');
    nElements = length(prependElements);
    nOut = nIn + nElements;
    outCell = cell(nOut,1);
    outCell(1:nElements) = prependElements(:);
    outCell((nElements+1):end) = inCell(:);
else
    nOut = nIn + 1;
    outCell = cell(nOut,1);
    outCell{1} = prependElements;
    outCell(2:end) = inCell(:);
end
end
    
%% curly
function out = curly(in,varargin)
%https://www.mathworks.com/matlabcentral/fileexchange/39735-functional-programming-constructs
%https://blogs.mathworks.com/loren/2013/01/24/introduction-to-functional-programming-with-anonymous-functions-part-2/
%map     = @(val, fcns) cellfun(@(f) f(val{:}), fcns);
%mapc    = @(val, fcns) cellfun(@(f) f(val{:}), fcns, 'UniformOutput', 0);
%iif     = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
%recur   = @(f, varargin) f(f, varargin{:});
%paren   = @(x, varargin) x(varargin{:});
%curly   = @(x, varargin) x{varargin{:}};
%
%NB: curly(aCellArray,':'); %is like aCellArray{:}

out = in{varargin{:}};
end

%% nanzscore
function zScoreData = nanzscore(data)
% column wise zscore with omitnan

zScoreData = (data - mean(data,1,'omitnan'))./std(data,0,1,'omitnan');
end

%% checkAllVariablesOI
function [figsOutlier,figsCorr,corrStruct] = checkAllVariablesOI(dataROIs,roiNames,regressorsTable,plotLevel)
% do basic check, outlier plots and correlations

%% variables
Sex           = regressorsTable.Sex;
Handedness    = regressorsTable.Handedness;
% SexHandedness = Sex.*Handedness;
Age           = regressorsTable.Age; %demean(regressorsTable.Age); %
dt_RSfMRI     = regressorsTable.dt_RSfMRI; %demean(regressorsTable.dt_RSfMRI); %

% ReidsL        = regressorsTable.ReidsL;
% ReidsR        = regressorsTable.ReidsR;
% ReidsAv       = regressorsTable.ReidsAv; %(ReidsL+ReidsR)./2;
% 
% angleMVS      = regressorsTable.angleMVS;
% angleMVSL     = regressorsTable.angleMVSL;
% angleMVSR     = regressorsTable.angleMVSR;

%transformed angles
pReidsL       = regressorsTable.pReidsL;
pReidsR       = regressorsTable.pReidsR;
pReidsAv      = regressorsTable.pReidsAv; 

pMVS          = regressorsTable.pMVS; 
pMVSL         = regressorsTable.pMVSL; 
pMVSR         = regressorsTable.pMVSR; 

%% some raw covariates plots: first find outliers then plot data and mark outliers
disp('checking for outliers in the "Angle" covariates...');
covarNames = {'pReidsL';'pReidsR';'pReidsAv';'pMVSL';'pMVSR';'pMVS'}; 
covarsOI   = [ pReidsL,  pReidsR,  pReidsAv,  pMVSL,  pMVSR,  pMVS]; 

allOutliers = zeros(size(covarsOI));
for indVar = 1:length(covarNames)
    currVarName = covarNames{indVar};
    currVar = covarsOI(:,indVar);
    combOutliers = [isoutlier(currVar,'mean'),isoutlier(currVar,'median'),isoutlier(currVar,'quartiles'),isoutlier(currVar,'grubbs'),isoutlier(currVar,'gesd')];
    sumCombOutliers = sum(combOutliers,2);
    allOutliers(:,indVar) = sumCombOutliers>3;
    if(any(allOutliers(:,indVar)))
        disp(['There are outliers in "',currVarName,'"! Will remove those for the regression. (robustfit can deal with nans.)']);
        %keyboard;
    else
        disp(['It seems "',currVarName,'" does NOT have outliers.']);
    end
end

figsOutlier = cell(3,1); %init
if(strcmp(plotLevel,'all')||strcmp(plotLevel,'sparse'))
figsOutlier{1} = figure(10); clf;
[S,AX,BigAx,H,HAx] = plotmatrix(covarsOI); title('pReids & pMVS'); hold('on')
for indP = 1:size(AX,1)
    AX(indP,1).YLabel.String = covarNames{indP};
    AX(end,indP).XLabel.String = covarNames{indP};
end
figsOutlier{2}.S_pmat  = S;
figsOutlier{2}.AX      = AX;
figsOutlier{2}.BigAx   = BigAx;
figsOutlier{2}.HistObj = H;
figsOutlier{2}.HistAx  = HAx;

%% now plot outliers again in red and normal data in blue
figsOutlier{3} = figure(11); clf;
for indCV2 = 1:length(covarNames)
    for indCV1 = 1:length(covarNames)
        subplot(length(covarNames),length(covarNames),((indCV2-1)*length(covarNames))+indCV1);
        currXLim = AX(indCV2,indCV1).XLim;
        currYLim = AX(indCV2,indCV1).YLim;
        if(indCV1==indCV2)
            histogram(covarsOI(:,indCV2),H(indCV2).BinEdges);
            xlim(currXLim);
            ylim([min(H(indCV2).Values); max(H(indCV2).Values)]);
        else
            currOutliers = (allOutliers(:,indCV1) + allOutliers(:,indCV2))~=0;
            plot(covarsOI(~currOutliers,indCV1),covarsOI(~currOutliers,indCV2),'bo','MarkerFaceColor','b'); hold('on');
            if(any(currOutliers))
                plot(covarsOI( currOutliers,indCV1),covarsOI( currOutliers,indCV2),'rx');
            end            
            xlim(currXLim);
            ylim(currYLim);
        end
        if(indCV1>1&&indCV2<length(covarNames))
            set(gca,'xtick',[],'ytick',[])
        end
        if(indCV1==1)
            ylabel(covarNames{indCV2});
            if(indCV2~=length(covarNames))
                set(gca,'xtick',[])
            end
        end
        if(indCV2==length(covarNames))
            xlabel(covarNames{indCV1})
            if(indCV1~=1)
                set(gca,'ytick',[])
            end
        end
        
    end
end
sgtitle('pReids & pMVS with outliers marked as red crosses');
end

%% check "raw" correlations for all covariates (in this new analysis)
disp('plotting correlations (just as a reference/additional information, will use robust regression later)...');
covarsMat  = [pMVSL,pMVSR,pMVS,Age,dt_RSfMRI,pReidsL,pReidsR,pReidsAv]; 
covarNames = {'pMVSL';'pMVSR';'pMVS';'Age';'dt_RSfMRI';'pReidsL';'pReidsR';'pReidsAv'}; 
[figsCorr,corrStruct] = checkPlotCorrs(roiNames,dataROIs,regexprep(covarNames,'_','-'),covarsMat,1:2);

end