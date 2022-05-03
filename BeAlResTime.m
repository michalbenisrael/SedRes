function [corrStart,corrEnd,expoAgeBeStart,expoAgeBeEnd] = BeAlResTime...
    (Be, Al,max_nBe_inhert,maxNumOfSteps,maxTimeOfRun,...
    burialLayerDepth,burialTimeMax,HdistType,TdistType,gamma_sp,gamma_mu,rho,LAl,...
    LBe,P_spBe,P_muBe,P_spAl,P_muAl,isplot)

% BeAlResTime calculates residence time of sediments in rivers using 10Be and 26Al data from large rivers and calculating how long it would take for these values to be reached during transport in the river. 
corrStart = NaN;
corrEnd = NaN;
expoAgeBeStart = NaN;
expoAgeBeEnd = NaN;

nBe_inhert = rand(1)*max_nBe_inhert; %this is the max number of Be atoms that we can deal with.
% the number of "max_nBe_inhert" should be derived from the minimum Be in
% the samples, and reflects the erosion rate in the region.
nAl_inhert = nBe_inhert*6.75;%(30.26/4.96); % Similarly, Al is the same as Be, to the factor of 6.

nAl0 = nAl_inhert;
nBe0 = nBe_inhert;

% depth
if HdistType==1
    h = exprnd(burialLayerDepth,[1,maxNumOfSteps]); %exp. depth dist.
elseif HdistType ==2
    h = rand([1,maxNumOfSteps]).*burialLayerDepth*2;
end

% dist.
if TdistType==1
    t = round(exprnd(burialTimeMax,[1,maxNumOfSteps])); %exp. age dist.
elseif TdistType==2
    t = rand([1,maxNumOfSteps]).*burialTimeMax;
end

tSum = cumsum(t);
close all

n = find(tSum>maxTimeOfRun,1); %finding how many model runs we need
if ~isempty(n)
    maxNumOfSteps = n; %if we can limit - we limit the number of steps
end

t = t(1:maxNumOfSteps);
tSum = tSum(1:maxNumOfSteps);
h = h(1:maxNumOfSteps);

nAl = nan(1,maxNumOfSteps); nAl(1) = nAl0;
nBe = nan(1,maxNumOfSteps); nBe(1) = nBe0;

i = 2;
while (i < maxNumOfSteps) && (nBe(i-1) < (Be(1)+Be(2)))
    nAl(i) = nAl(i-1)*exp(-LAl*t(i)) + (( P_spAl*exp(-(rho/gamma_sp)*h(i)) + P_muAl*exp(-(rho/gamma_mu)*h(i)) )/LAl)*(1-exp(-LAl*t(i)));
    nBe(i) = nBe(i-1)*exp(-LBe*t(i)) + (( P_spBe*exp(-(rho/gamma_sp)*h(i)) + P_muBe*exp(-(rho/gamma_mu)*h(i)) )/LBe)*(1-exp(-LBe*t(i)));
    i=i+1;
end

AlCorr = ismembertol(nAl,Al(1),Al(2)/max(abs([nAl(:);Al(1)])));
BeCorr = ismembertol(nBe,Be(1),Be(2)/max(abs([nBe(:);Be(1)])));
iso2 = find(BeCorr+AlCorr==2);
if ~isempty(iso2)
    corrStart = tSum(iso2(1));
    corrEnd = tSum(iso2(end));
    expoAgeBeStart = (nBe(iso2(1))-nBe0)/P_spBe;
    expoAgeBeEnd = (nBe(iso2(end))-nBe0)/P_spBe;
else
    return
end
