%% Phase Flip Analysis Pipeline
clear 
clc
subj=5;

% Load
cd 'Y:\el-Christina\PhaseFlip\PF_Pilot\Raw Data\Raw Behaviour'
loadfilename=sprintf('Pilot_PhaseFlip_Subj%i.mat',subj);
load(loadfilename)

%% JND Task
figure;
JND=mean(subresults.JND_UD.x(end-10:end));
plot(1:length(subresults.JND_results.Comparison), subresults.JND_results.Comparison, '-o', 'LineWidth', 2);
xlim([1 length(subresults.JND_results.Comparison)]);
ylim([min(subresults.JND_results.Comparison)-1 max(subresults.JND_results.Comparison)+1]);
title(sprintf('JND Results (JND: %.3f)',JND));
xlabel('Trial Number');
ylabel('Comparison Duration');
yline(JND);

%% Gabor intensity progression for each condition separately
gabor_cond1=subresults.data{subresults.data{:,'Condition'}==1,'Gabor Strength'};
gabor_cond2=subresults.data{subresults.data{:,'Condition'}==2,'Gabor Strength'};

figure;
plot(1:length(gabor_cond1),gabor_cond1, '-o', 'LineWidth', 2);
hold on
plot(1:length(gabor_cond1),gabor_cond2, '-o', 'LineWidth', 2);
title('Gabor Strength');
xlabel('Trial Number');
ylabel('Gabor Intensity');
box('off')

%% Performance valid vs invalid
data_cond1=subresults.data(subresults.data{:,'Condition'}==1,:);
data_cond2=subresults.data(subresults.data{:,'Condition'}==2,:);

% Get correct/incorrect response data for different timings and conditions
short_cond1=data_cond1{data_cond1{:,'Target Interval'}==0.5,'Correct/Incorrect'};
long_cond1=data_cond1{data_cond1{:,'Target Interval'}==1.2,'Correct/Incorrect'};
valid_cond1=data_cond1{data_cond1{:,'Target Interval'}==0.7,'Correct/Incorrect'};

short_cond2=data_cond2{data_cond2{:,'Target Interval'}==0.55,'Correct/Incorrect'};
long_cond2=data_cond2{data_cond2{:,'Target Interval'}==1.25,'Correct/Incorrect'};
valid_cond2=data_cond2{data_cond2{:,'Target Interval'}==0.75,'Correct/Incorrect'};

% Get means
short_cond1=mean(short_cond1);
long_cond1=mean(long_cond1);
valid_cond1=mean(valid_cond1);

short_cond2=mean(short_cond2);
long_cond2=mean(long_cond2);
valid_cond2=mean(valid_cond2);

% Plot by validity
figure;
bar([short_cond1 short_cond2; long_cond1 long_cond2; valid_cond1 valid_cond2])
xticklabels(["Short","Long","Valid"])
legend(["800","850"])
title("1 Performance - Valid vs. Invalid")
ylabel('Percent Correct')
box('off')

% Plot by conditions
figure;
bar([short_cond1 long_cond1  valid_cond1; short_cond2 long_cond2 valid_cond2])
xticklabels(["800","850"])
legend(["Short","Long","Valid"])
title("2 Performance - Valid vs. Invalid")
ylabel('Percent Correct')
box('off')

%%