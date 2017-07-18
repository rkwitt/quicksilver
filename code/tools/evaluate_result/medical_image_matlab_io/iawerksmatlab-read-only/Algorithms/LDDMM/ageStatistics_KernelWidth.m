close all; clear all


% where to truncate gaussians
weightCutoff = 0.04;

% sigma for gaussians
sigmaValues = [1 2 4:4:32];
%sigmaValues = 4:0.5:10;

% ages where densities will be computed
regressionAges = -35:1:125;

% age values for atlas centers
atlasCenterAges = 20:1:80;

% where (on the y axis) to align the age scatter plot
yAlignmentAgeScatterplot = 0.05;

% ID numbers of the male patients; ordering correspoinds to ages below
idsmale = [
     6
    34
     2
    35
    13
    21
    26
    12
    19
    55
    37
    79
    57
    59
    80
    88
    45
     7
    53
    10
    39
    63
    93
    56
    89
    86
    58
    44
    25
    94
    81
    68
    91
    30
    95
    92
    46
   100
    99
    41
    97
    29
   105
];

% ages of the male patients; ordering corresponds to ids above
agesmale = [
    22
    23
    25
    25
    26
    26
    26
    28
    28
    30
    31
    31
    33
    33
    34
    34
    35
    36
    37
    41
    41
    43
    43
    44
    44
    46
    47
    48
    49
    50
    51
    56
    56
    57
    59
    60
    61
    61
    63
    64
    68
    72
    79
];

% ID numbers of the female patients; ordering correspoinds to ages below
idsfemale = [
%    003   % found to not exist
    082 
    009 
    033 
    005 
    040 
    020 
    004 
    011 
    050 
    047 
    038 
    071 
    008 
    062 
    054 
    064 
    042 
    022 
    017 
    023 
    074 
    018 
    065 
    052 
    027 
    060 
    070 
    031 
    072 
    085 
    051 
    024 
    061 
    067 
    087 
    090 
    075 
    073 
    036 
    049 
    069 
    028 
    032 
    066 
    096
];

% ages of the female patients; ordering corresponds to ids above
agesfemale = [
%    20    % found to not exist
    22 
    23 
    23 
    24 
    25 
    26 
    27 
    27 
    28 
    31 
    32 
    32 
    33 
    33 
    34 
    35 
    37 
    39 
    40 
    40 
    41 
    43 
    43 
    44 
    48 
    48 
    48 
    50 
    52 
    52 
    55 
    57 
    57 
    57 
    57 
    57 
    61 
    64 
    65 
    65 
    65 
    66 
    66 
    66 
    68
];

% create a list of all patients
idAll = [idsmale; idsfemale];
agesAll = [agesmale; agesfemale];

[agesAll, sortIndex] = sort(agesAll);
idAll = idAll(sortIndex);

% create a list of shuffled ages
perm = randperm(length(idAll));
agesAllRandom = agesAll(perm);

% jitterplot of patients
%plot(agesAll,...
%       yAlignmentAgeScatterplot * ones(1,length(agesAll)) + ...
%       0.03*(rand(1,length(agesAll))-0.5),...
%       'ok','LineWidth',2);
%hold on;

kernelDensities = zeros(length(sigmaValues), length(regressionAges));
currentWeights  = zeros(length(agesAll),1);

% kernel density estimate for each sigma
for sigmaIndex = 1:length(sigmaValues)
  for centerIndex = 1:length(regressionAges)
    for observationIndex = 1:length(agesAll)
      currentWeights(observationIndex) = ...
        regkernel(agesAll(observationIndex), ...
        regressionAges(centerIndex), ...
        sigmaValues(sigmaIndex));
    end
    kernelDensities(sigmaIndex,centerIndex) = mean(currentWeights);
  end
end

%plot(regressionAges,kernelDensities);
%legendStrs = {'Observations'};
%for s = 1:length(sigmaValues)
%  legendStrs{s+1} = sprintf('sigma=%s',num2str(sigmaValues(s)));
%end
%legend(legendStrs);
%title({'Kernel Density Estimation For Varying Sigma Values','Population Size = 88'});
%xlabel('Age');
%ylabel('Density');

%
% make a plot of an example kernel 
%
%f = figure;
% jitterplot of patients
%plot(agesAll,...
%       yAlignmentAgeScatterplot*3 * ones(1,length(agesAll)) + ...
%       0.03*(rand(1,length(agesAll))-0.5),...
%       '+k','LineWidth',2);
%hold on;

centerAge = 44;
ctsAges = -40:0.1:130;
exampleDensities = zeros(length(sigmaValues),length(ctsAges));
for s = 1:length(sigmaValues)
  for i = 1:length(ctsAges)
    exampleDensities(s,i) = regkernel(ctsAges(i), ...
      centerAge, sigmaValues(s));
  end
end
%plot(ctsAges(1,:),exampleDensities(1,:), '-k','LineWidth',2);
%plot(ctsAges(1,:),exampleDensities(2,:), ':k','LineWidth',2);
%plot(ctsAges(1,:),exampleDensities(3,:), '--k','LineWidth',2);
%title('Illustration of Regression Kernel Weights');
%xlabel('Patient Age');
%ylabel('Weight');
%legend(legendStrs);
%xlim([centerAge-25 centerAge+25]);
%ylim([0 0.25]);
%ca = gca;
%boldifyPlot(f,ca);

%
% now compute actual weights
%
outfid = fopen('Weights_25April07.txt','w');
weights = zeros(length(sigmaValues),length(atlasCenterAges),length(agesAll));
atlasSizeTable = zeros(length(sigmaValues),length(atlasCenterAges));
for sigmaIndex = 1:length(sigmaValues)
  for centerIndex = 1:length(atlasCenterAges)
    for observationIndex = 1:length(agesAll)
      weights(sigmaIndex,centerIndex,observationIndex) = ...
        regkernel(agesAll(observationIndex), ...
        atlasCenterAges(centerIndex),...
        sigmaValues(sigmaIndex));
    end
    
    % print how many images will be need to create each atlas
    w = weights(sigmaIndex,centerIndex,:);
    w = w(:);

    m1Cutoff = regkernel(1*sigmaValues(sigmaIndex), 0, sigmaValues(sigmaIndex));
    m2Cutoff = regkernel(2*sigmaValues(sigmaIndex), 0, sigmaValues(sigmaIndex));
    m3Cutoff = regkernel(3*sigmaValues(sigmaIndex), 0, sigmaValues(sigmaIndex));
      
    m1 = sum(w>m1Cutoff);
    m2 = sum(w>m2Cutoff);
    m3 = sum(w>m3Cutoff);
    atlasSizeTable(sigmaIndex,centerIndex) = m3;
    fprintf(outfid, 'Std. Dev: %d, Atlas Age: %d, M1: %d, M2: %d, M3: %d\n',...
      sigmaValues(sigmaIndex), atlasCenterAges(centerIndex), m1, m2, m3);

    % three standard deviations cuttoff
    atlasIDs      = idAll(w>m3Cutoff);
    kernelWeights = w(w>m3Cutoff);
    atlasAges     = agesAll(w>m3Cutoff);

    ageDifferences = atlasAges - atlasCenterAges(centerIndex);
    
    p0Weights     = kernelWeights./sum(kernelWeights(:));
    
    s0 = 1/m3 * sum(kernelWeights .* ageDifferences.^0);
    s1 = 1/m3 * sum(kernelWeights .* ageDifferences.^1);
    s2 = 1/m3 * sum(kernelWeights .* ageDifferences.^2);
    
    tmpNumerator = s2 - s1*ageDifferences;
    tmpDenomonator = s2*s0 - s1^2;
    
    p1Weights = 1/m3 * kernelWeights .* tmpNumerator / tmpDenomonator;
    
    fprintf(outfid,'Sum P0: %g, Sum P1: %g\n', sum(p0Weights), sum(p1Weights));
    fprintf(outfid,'%03d  %03d %g %g %g\n', [atlasIDs atlasAges kernelWeights p0Weights p1Weights]');
  end
end
fclose(outfid);

%
% comute weights for randomized data
%
weightsRandom = zeros(length(sigmaValues),length(atlasCenterAges),length(agesAll));
for sigmaIndex = 1:length(sigmaValues)
  for centerIndex = 1:length(atlasCenterAges)
    for observationIndex = 1:length(agesAll)
      weightsRandom(sigmaIndex,centerIndex,observationIndex) = ...
        regkernel(agesAllRandom(observationIndex), ...
        atlasCenterAges(centerIndex),...
        sigmaValues(sigmaIndex));
    end
  end
end

atlasSizeTable