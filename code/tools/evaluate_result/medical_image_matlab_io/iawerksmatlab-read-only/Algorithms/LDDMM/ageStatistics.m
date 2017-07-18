close all; clear all

% hack to work around plotting but in octave
plotMale = true;
plotFemale = true;

% where (on the y axis) to align the age scatter plot
yAlignmentAgeScatterplot = 0.135;
% where (on the y axis) to align the age scatterplot for each kernel
yAlignmentKernelAgeScatterplot = 0.2:0.01:0.3;

% where to truncate gaussians
weightCutoff = 0.04;

% sigma for gaussians
sigmaValues = 5.0;

% age values for atlas centers
atlasCenterAges = 35:4:55;

figureMale = figure;
figureFemale = figure;

outfid = fopen('Weights_Sigma5.txt','w');

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

%
% create a scatterplot of the patients by age, align the scatterplots
% at y=yAlignmentAgeScatterplot
%
if plotMale
  figure(figureMale); hold on;
  plot(agesmale,...
       yAlignmentAgeScatterplot * ones(1,length(agesmale)) + ...
       0.03*(rand(1,length(agesmale))-0.5),...
       'ok','LineWidth',2);
end

if plotFemale
  figure(figureFemale); hold on;
  plot(agesfemale,...
       yAlignmentAgeScatterplot * ones(1,length(agesfemale)) + ...
       0.03*(rand(1,length(agesfemale))-0.5),...
       'ok','LineWidth',2);
end

%
% compute density at each age given the current kernel
%
densityMale = zeros(length(sigmaValues),length(atlasCenterAges));
densityFemale = zeros(length(sigmaValues),length(atlasCenterAges));

%
% store the error between theoretical and actual kernels for each
% sigma
%
errorArrayMale   = zeros(length(sigmaValues),1);
errorArrayFemale = zeros(length(sigmaValues),1);

tages = 20:0.1:70;
theoreticalKernel = zeros(length(tages));
for i = 1:length(tages)
  theoreticalKernel(i) = regkernel(tages(i),45, sigmaValues(1));
end

for s=1:length(sigmaValues)
  currentSigma = sigmaValues(s);

  for j=1:length(atlasCenterAges)
    atlasAge = atlasCenterAges(j);

    %
    % males
    %

    % compute theoretical weights
    weights = zeros(1,length(agesmale));
    for i = 1:length(agesmale)
      weights(i) = regkernel(agesmale(i),atlasAge,currentSigma);
      densityMale(s,j) = densityMale(s,j) + weights(i);
    end

    % normalize weights
    weights = weights / sum(weights(:));
    % plot original kernel weights
    if plotMale
      figure(figureMale);
      plot(agesmale, weights, '-k','LineWidth',2.5);
      title({'Kernel Weights: Ages 35-55','Males'});
      xlabel('Age (years)');
      ylabel('Weight');
      ylim([0 0.2]); 
      xlim([15 75]);
      %legend('Image',');
    end
    weightsOrig = weights;

    % truncate kernel and renormalize
    weights(weights < weightCutoff) = 0;
    weights = weights / sum(weights(:));
    % plot truncated kerenl weights
    if plotMale
%      plot(agesmale(weights>0), weights(weights>0), '-');
    end

    % compute error of truncated kernel
    fitError = sum((weightsOrig - weights).^2);
    errorArrayMale(s) = errorArrayMale(s) + fitError;

    % scatterplot of support for this kernel
    if false & plotMale
      plot(agesmale(weights>0),...
           yAlignmentKernelAgeScatterplot(mod(j,length(yAlignmentKernelAgeScatterplot))+1) * ...
           ones(1,length(agesmale(weights>0))) + ...
           (rand(1,length(agesmale(weights>0)))-0.5)*0.0,'-');
    end

    % print id's and corresponding non-zero weights 
%    fprintf(outfid, 'Male, Age=%d, M=%d, Sigma=%g, Sum Weights=%g\n',...
%            atlasAge, sum(weights > 0),currentSigma, sum(weights));
%    fprintf(outfid, '%03d  %g\n', [idsmale(weights > 0)'; weights(weights > 0)](:));
%    fprintf(outfid, '\n');
    %[idsmale(weights > 0) weights(weights > 0)']

    %
    % Females
    %

    % compute theoretical weights
    weights = zeros(1,length(agesfemale));
    for i = 1:length(agesfemale)
      weights(i) = regkernel(agesfemale(i),atlasAge,currentSigma);
      densityFemale(s,j) = densityFemale(s,j) + weights(i);
    end

    % normalize weights
    weights = weights / sum(weights(:));
    % plot original kernel weights
    if plotFemale
      figure(figureFemale);
      plot(agesfemale, weights, '-k', 'LineWidth',2.5);
      title({'Kernel Weights: Ages 35-55','Females'});
      xlabel('Age (years)');
      ylabel('Weight');   
      ylim([0 0.2]);
      xlim([15 75]);

      %legend('Image','Emperical Weights', 'Theoretical Weights');
    end
    weightsOrig = weights;

    % truncate kernel and renormalize
    weights(weights < weightCutoff) = 0;
    weights = weights / sum(weights(:));
    % plot truncated kerenl weights
    if plotFemale
%      plot(agesfemale(weights>0), weights(weights>0), '-');
    end

    % compute error of truncated kernel
    fitError = sum((weightsOrig - weights).^2);
    errorArrayFemale(s) = errorArrayFemale(s) + fitError;

    % scatterplot of support for this kernel
    if false & plotFemale
      plot(agesfemale(weights>0),...
           yAlignmentKernelAgeScatterplot(mod(j,length(yAlignmentKernelAgeScatterplot))+1) * ...
           ones(1,length(agesfemale(weights>0))) + ...
           (rand(1,length(agesfemale(weights>0)))-0.5)*0.0,'-');
    end

    % print id's and corresponding non-zero weights 
%    fprintf(outfid, 'Female, Age=%d, M=%d, Sigma=%g, Sum Weights=%g\n',...
%            atlasAge, sum(weights > 0),currentSigma, sum(weights));
%    fprintf(outfid, '%03d  %g\n', [idsfemale(weights > 0)'; weights(weights > 0)](:));
%    fprintf(outfid, '\n');
    %[idsfemale(weights > 0) weights(weights > 0)']
  end
end
%figure(figureMale);
%plot(tages,theoreticalKernel,'-r','LineWidth',2.5);
%figure(figureFemale);
%plot(tages,theoreticalKernel,'-r','LineWidth',2.5);

%
% plot the densities
%
if false & plotMale
  figure(figureMale);
  plot(atlasCenterAges, yAlignmentAgeScatterplot*densityMale', '-b');
  title(sprintf('Male Kernel Weights (Truncated & Normalized) sig=%g', sigmaValues(1)))
  xlabel('Age');
  ylabel('Kernel Weight');
end

if false & plotFemale
  figure(figureFemale);
  plot(atlasCenterAges,yAlignmentAgeScatterplot*densityFemale','-b');
  title(sprintf('Female Kernel Weights (Truncated & Normalized) sig=%g', sigmaValues(1)))
  xlabel('Age');
  ylabel('Kernel Weight');
end

%figureError = figure;
%plot(sigmaValues', [errorArrayMale errorArrayFemale]);
%legend('Male', 'Female')
%title('SSE Theoretical-Truncated Kernel Fit');
%xlabel('sigma');
%ylabel('sse');
boldifyPlot(figureMale, gca);
boldifyPlot(figureFemale, gca);
fclose(outfid);