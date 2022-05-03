clc, clear, close all
% This is the father function that where variables and constants a defined for the main function “BeAlResTime”
% All nuclide concentration and their uncertainties in atoms/g

% Relevant data for the different rivers
% %Po River
% River='Po River';
% Be(1,:)=[3.30, 3.66, 2.98].*1e4;
% Be(2,:)=[0.16, 0.19, 0.25].*1e4;
% Al(1,:)=[20.2, 23.6, 17.4].*1e4;
% Al(2,:)=[2.1, 2.7, 2.1].*1e4;
% lat = 45; elv = 50; %Averaged elevation and latitude for Po River
% samplesName={'P1' 'P3-1' 'P6-2'};  %Po
%  burialdepth = 20;
% burialtime = 20;

% %Branco River (Amazon tributary)
% River='Branco (Amazon tributary)';
% Be(1,:)=[35.10 31.99 46.91 33.87 32.20].*1e4;
% Be(2,:)=[1.56 1.55 3.95 2.22 2.10].*1e4;
% Al(1,:)=[216.03 185.64 183.26 189.78 169.62].*1e4;
% Al(2,:)=[15.15 16.74 17.16 12.01 12.14].*1e4;
% lat = 1; elv = 50; %Averaged elevation and latitude for Branco River
% samplesName={'Br1a' 'Br4c' 'Br5b' 'Br5c' 'Br8b2'};  %Branco
%  burialdepth = 2;
% burialtime = 20;

% %Amazon lowlands
% River='Amazon lowlands';
% Be(1,:)=[7.60 8.01 8.77 13.31 6.2 8.32].*1e4;
% Be(2,:)=[0.33 0.69 0.56 0.47 0.36 0.4].*1e4;
% Al(1,:)=[32.81 48.09 49.36 92.86 36.49 34.6].*1e4;
% Al(2,:)=[8.74 4.66 16.18 8.26 4.14 10.27].*1e4;
% lat = 1; elv = 5; %Averaged elevation and latitude for Amazon lowlands
% samplesName={'Man1.1b' 'Man1.1c-2' 'Ir0.4b' 'Ir0.4c' 'Par0.9a' 'Ama-b'};  
%  burialdepth = 6;
% burialtime = 20;

% %Colorado
% River='Lower Colorado';
% Be(1,:)=[1.17 1.71 9.75 2.85 9.89 1.33].*1e5;
% Be(2,:)=[0.983 0.887 0.945 1.48 0.959 1.27].*1e4;
% Al(1,:)=[6.67 9.06 5.79 15.1 6.51 7.55].*1e5; 
% Al(2,:)=[1.62 2.20 1.41 3.68 1.58 1.83].*1e5; 
% lat = 33; elv = 100; %Averaged elevation and latitude for lower Colorado
% samplesName={'WB' 'TC' 'ND' 'PD' 'ER' 'YM'}; %Colorado
%  burialdepth = 10; 
% burialtime = 20;


gamma_sp = 1600; % kg/m2
gamma_mu = 15000; % kg/m2
rho = 2200; % kg/m3
PBe=4.96; FspBe=0.978;% Production rate at SLHL
PAl=30.26; FspAl=0.974;
LAl=9.79E-07; LBe=4.997E-07; %Decay const for Al-Be 1/yr
%LAl=9.83E-07; LBe=4.62E-07; %Decay const for Al-Be 1/yr
%production is calculated based on procedures presented on Stone (2000)
[P_spAl,P_muAl]=Stone2000(elv,lat,PAl,FspAl);
[P_spBe,P_muBe]=Stone2000(elv,lat,PBe,FspBe);

max_nBe_inhert = min(Be(1,:)-Be(2,:)); % maximum number of inherited Be atoms

maxNumOfSteps = 1e6; %[maximum # of model steps]
maxTimeOfRun = 1e6;
%burial depth and time change depending on data from each river
burialLayerDepth = burialdepth; %[m] changes per river
burialTimeMax = burialtime; %[yr] constant for rivers


% See paper text for reasoning behind H and t distribution types
HdistType = 2; % 1 - exponential distribution (parameter is the mean of the dist.); 2 - Uniform distribution (parameter is the max).
TdistType = 1; % likewise for the times

if HdistType==1
    if TdistType==1
        plotsDir = strcat(plotsDir,'expT_expH\');
    elseif TdistType==2
        plotsDir = strcat(plotsDir,'uniT_expH\');
    end
elseif HdistType==2
     if TdistType==1
        plotsDir = strcat(plotsDir,'expT_uniH\');
    elseif TdistType==2
        plotsDir = strcat(plotsDir,'uniT_uniH\');
     end
end


for j = 1:length(Be)
    sample = j;
    
    for i = 1:1000
        [corrStart(j,i),corrEnd(j,i),expoAgeBeStart(j,i),expoAgeBeEnd(j,i)] = BeAlResTime...
        ([Be(1,sample), Be(2,sample)], [Al(1,sample),Al(2,sample)],max_nBe_inhert,maxNumOfSteps,maxTimeOfRun,...
        burialLayerDepth,burialTimeMax,HdistType,TdistType,gamma_sp,gamma_mu,rho,LAl,...
        LBe,P_spBe,P_muBe,P_spAl,P_muAl,0);
   
    end
    
end
%%

ResMed=median(corrAll,'omitnan');
ResSTD=std(corrAll,'omitnan');
ResidenceTime(1)=max(ResMed,[],'omitnan');
ResidenceTime(2)=min(ResMed,[],'omitnan');
