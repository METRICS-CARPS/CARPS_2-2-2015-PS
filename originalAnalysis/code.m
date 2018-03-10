function code
% this code will approximate results from Sullivan, N. J., Hutcherson, C. A., 
% Harris, A., Rangel, A. Dietary self-control is related to the speed with which 
% health and taste attributes are processed. Psychological Science. Forthcoming.
%
% written 10/2014 by Nikki Sullivan

%% load and format data

% x, y dimensions of computer screen:
screenDims = [1680 1050];
cursorStartPixels = [840 840];
subjects=1:28;

% get mouse data
% columns in the csv file:
% [subject | time point | x-position (pixels) | y-position (pixels) | 
% liking rating (left) | liking rating (right) | 
% taste rating (left) | taste rating (right) | 
% health rating (left) | health rating (right) |choice (1==left) | RT ]
data=csvread('sullivanEtAl2014ForPub_mouse.csv');
nTPs = NaN(280,28);
xPos = NaN(280,8500,28);
yPos = NaN(280,8500,28);
like = NaN(280,2,28);
taste = NaN(280,2,28);
health = NaN(280,2,28);
choice = NaN(280,28);
RT = NaN(280,28);
healthChosen = NaN(280,28);
healthUnchosen = NaN(280,28);
unchosen=[2 1];
for subj = 1:28
    % find this subject
     thisSubjData = data(data(:,1)==subj,:);
    
    % find each trial
    trialStarts = find(thisSubjData(:,2)==1);
    trialEnds = [trialStarts(2:end)-1; length(thisSubjData)]; 
    
    % separate data by trial
    for trial = 1:length(trialStarts)
        nTPs(trial,subj) = (1+trialEnds(trial)) - trialStarts(trial);
        
        xPos(trial,1:nTPs(trial,subj),subj) = thisSubjData(trialStarts(trial):trialEnds(trial),3);
        yPos(trial,1:nTPs(trial,subj),subj) = thisSubjData(trialStarts(trial):trialEnds(trial),4);
        
        like(trial,:,subj) = thisSubjData(trialStarts(trial),5:6);
        taste(trial,:,subj) = thisSubjData(trialStarts(trial),7:8);
        health(trial,:,subj) = thisSubjData(trialStarts(trial),9:10);
        choice(trial,subj) = thisSubjData(trialStarts(trial),11);
        RT(trial,subj) = thisSubjData(trialStarts(trial),12);
        
        healthChosen(trial,subj) = health(trial,choice(trial,subj),subj);
        healthUnchosen(trial,subj) = health(trial,unchosen(choice(trial,subj)),subj);
    end
    
end

% get keyboard data
% columns:
% [subject |  
% liking rating (left) | liking rating (right) | 
% taste rating (left) | taste rating (right) | 
% health rating (left) | health rating (right) | choice | RT ]
data=csvread('sullivanEtAl2014ForPub_key.csv');
key_like = NaN(80,2,28);
key_taste = NaN(80,2,28);
key_health = NaN(80,2,28);
key_healthChosen = NaN(80,28);
key_healthUnchosen = NaN(80,28);
key_choice = NaN(80,28);
key_RT = NaN(80,28);
for subj = 1:28
    % find this subject
    thisSubjData = data(data(:,1)==subj,:);
    nTrials = length(thisSubjData);
    
    key_like(1:nTrials,1:2,subj) = thisSubjData(:,2:3);
    key_taste(1:nTrials,1:2,subj) = thisSubjData(:,4:5);
    key_health(1:nTrials,1:2,subj) = thisSubjData(:,6:7);
    key_choice(1:nTrials,subj) = thisSubjData(:,8);
    key_RT(1:nTrials,subj) = thisSubjData(:,9);
    
    for trial = 1:nTrials
        key_healthChosen(trial,subj) = key_health(trial,key_choice(trial,subj),subj);
        key_healthUnchosen(trial,subj) = key_health(trial,unchosen(key_choice(trial,subj)),subj);
    end
end


%% transform mouse data


spat = cell(28,1);
spatDirect = cell(28,1);
temporal = cell(28,1);
temporalDirect = cell(28,1);
spatInterp = cell(28,1);
angleDirect = cell(28,1);
nTrials = NaN(28,1);

for subj = 1:28 % subject
    spat{subj} = NaN(8500,2,280);
    spatDirect{subj} = NaN(8500,2,280);
    temporal{subj} = NaN(101,2,280);
    temporalDirect{subj} = NaN(101,2,280);
    spatInterp{subj} = NaN(101,2,280);
    angleDirect{subj} = NaN(101,280);
    for trial = find(~isnan(choice(:,subj)))' % trial
        
        % spatial remapping: rescale from pixels to standard coord. space
        % convert all clicks to right-hand side - end is all (1,1)
        spat{subj}(1:nTPs(trial,subj),1,trial) = ...
            (xPos(trial,1:nTPs(trial,subj),subj) - xPos(trial,1,subj)) ...
            / (xPos(trial,nTPs(trial,subj),subj) - xPos(trial,1,subj));
        spat{subj}(1:nTPs(trial,subj),2,trial) = ...
            (yPos(trial,1:nTPs(trial,subj),subj) - yPos(trial,1,subj)) ...
            / (yPos(trial,nTPs(trial,subj),subj) - yPos(trial,1,subj));
        % retain choice - end is (-1,1) or (1,1)
        spatDirect{subj}(:,1,trial) = spat{subj}(:,1,trial) * sign(choice(trial,subj)-2+.5);
        spatDirect{subj}(:,2,trial) = spat{subj}(:,2,trial);

        % linear temporal{subj} interpolation (101 timepoints)
        timeStamps = (1:nTPs(trial,subj))';
        timeSamples = linspace(timeStamps(1),timeStamps(end),101);
        % convert all clicks to right-hand side - end is all (1,1)
        temporal{subj}(:,1,trial) = ...
            interp1(timeStamps,... % timepoints
            spat{subj}(1:nTPs(trial,subj),1,trial),... % coordinate
            timeSamples,... % timepoints to use for interploating
            'linear');
        temporal{subj}(:,2,trial) = interp1(timeStamps, ...
            spat{subj}(1:nTPs(trial,subj),2,trial),timeSamples,'linear');
        % retain choice - end is (-1,1) or (1,1)
        temporalDirect{subj}(:,1,trial) = ...
            interp1(timeStamps,spatDirect{subj}(1:nTPs(trial,subj),1,trial), ...
            timeSamples,'linear');
        temporalDirect{subj}(:,2,trial) = ...
            interp1(timeStamps,spatDirect{subj}(1:nTPs(trial,subj),2,trial), ...
            timeSamples,'linear');
        
        % temporal interpolation of trajectories with pixel info (fig. 2C)
        spatInterp{subj}(:,1,trial) = ...        
            interp1(timeStamps,xPos(trial,1:nTPs(trial,subj),subj), ...
            timeSamples,'linear');
        spatInterp{subj}(:,2,trial) = ...
            interp1(timeStamps,yPos(trial,1:nTPs(trial,subj),subj), ...
            timeSamples,'linear');

        
        % angle
        angleDirect{subj}(:,trial) = ...
            abs(atand(temporalDirect{subj}(:,1,trial) ./ ...
            temporalDirect{subj}(:,2,trial))) .* ...
            sign(temporalDirect{subj}(:,1,trial));
        % optional: exclude points where mouse moving downward
        excludeInd = [false; temporalDirect{subj}(2:end,2,trial) < ...
            temporalDirect{subj}(1:end-1,2,trial)];
        angleDirect{subj}(excludeInd,trial) = NaN;
        
    end
    nTrials(subj) = find(~isnan(squeeze(temporal{subj}(1,1,:))),1,'last');
end


%% excluding mouse trials


exclude = false(280,28);
reason.yCross =  false(280,28);
reason.RT =  false(280,28);
reason.overshoot =  false(280,28);
reason.smallRT =  false(280,28);
for subj = 1:28 % subject
    
    for trial = 1:nTrials(subj) % trial
        
        % how many times does cursor cross y axis? exclude if > 3
        ycross=0;
        for t=2:nTPs(trial,subj)
            if (roundn(xPos(trial,t,subj),1) < cursorStartPixels(1) && roundn(xPos(trial,t-1,subj),1) >= cursorStartPixels(1)) ...
                || (roundn(xPos(trial,t,subj),1) > cursorStartPixels(1) && roundn(xPos(trial,t-1,subj),1) <= cursorStartPixels(1))
                ycross = ycross+1;
            end
        end
        if ycross > 3
            exclude(trial, subj) = true;
            reason.yCross(trial, subj) = true;
        end
        
   end
    
    % exclude b/c RT > 1 or 2 SD above mean
    % in calculating standard deviation, don't include trials where RT > 5 min b/c that's a mistake or question they're asking (i.e. they're not doing the task).
    ind = ~exclude(:,subj) & RT(:,subj) <= 5;
    rtData = RT(ind,subj);
    oneSD = nanstd(rtData);
    twoSD = oneSD*2;
    ave = nanmean(rtData);
    excl.oneSD = RT(:,subj) > ave+oneSD;
    excl.twoSD = RT(:,subj) > ave+twoSD;
    output.excludeTrialOneSDRT = oneSD;
    output.excludeTrialTwoSDRT = twoSD;
    exclude(excl.twoSD,subj) = true;
    reason.RT(excl.twoSD,subj) = true;

    % **optional**: exclude "overshoots"
    % (times where the cursor goes far past food box before clicking)
    for trial = 1:nTrials(subj) % trial
        if sum(abs(spatDirect{subj}(:,1,trial)) >= 1.3) > 0
            exclude(trial,subj) = true;
            reason.overshoot(trial,subj) =true;
        end
    end 
    
    % **optional**: exclude very small RTs b/c that's weird
    if RT(:,subj) <= .2
        exclude(:, subj) = true;
        reason.smallRT(:, subj) = true;
    end
    
    % now remove those trials:
    angleDirect{subj}(:,exclude(:,subj)) = NaN;
    
end



%% self-control success ratio (SCSR) & control trial reaction time

successratio = zeros(1,28);
for subj = 1:28
    
    % classification uses both mouse and key trials:
    allHealth = [health(:,:,subj); key_health(:,:,subj)];
    allTaste = [taste(:,:,subj); key_taste(:,:,subj)];
    allHealthChosen = [healthChosen(:,subj); key_healthChosen(:,subj)];
    allHealthUnchosen = [healthUnchosen(:,subj); key_healthUnchosen(:,subj)];
    
    % challenge trial: trials requiring self control 
    % false = control not required. true=control required
    challengeTrial = ...
        ((allHealth(:,1) > allHealth(:,2)) & ...
        (allTaste(:,1) < allTaste(:,2))) | ...
        ((allHealth(:,2) > allHealth(:,1)) & ...
        (allTaste(:,2) < allTaste(:,1)));
    
    % healthier food chosen = true, otherwise false
    successTrial = challengeTrial & allHealthChosen > allHealthUnchosen;
    
    % get reaction times in each trial type
    allRT = [RT; key_RT];
    rtSCReq(subj) = nanmean(allRT(challengeTrial,subj));
    rtSCNotReq(subj) = nanmean(allRT(~challengeTrial,subj));
    
    % ratio
    successratio(subj) = sum(successTrial)/sum(challengeTrial);    

end

controlInd{1} = successratio > nanmean(successratio);
controlInd{2} = successratio < nanmean(successratio);


%% time course regressions

minTrials = 30; % minimum trials req'd for regression (only really important in non-normalized time)

l.coef = NaN(101,28);
l.pval = NaN(101,28);
th.coef = NaN(101,2,28);
th.pval = NaN(101,2,28);
for subj = 1:28 % subject

    mouseVar = angleDirect{subj};
    
    % liking difference regression
    rateDiff = like(:,2,subj) - like(:,1,subj);
    for t = 1:101 % one regression per time point
        if sum(~isnan(mouseVar(t,:))) > minTrials
            stats = regstats(mouseVar(t,:), rateDiff, 'linear', {'tstat'});
            l.coef(t,subj) = stats.tstat.beta(2);
            l.pval(t,subj) = stats.tstat.pval(2);
        end
    end
    % st calc: get one-sided t-test
    test = l.pval(:,subj)./2 < .05;
    l.st(subj) = 1+find(test==0,1,'last');
    
    % taste health difference regression
    rateDiff = [taste(:,2,subj) - taste(:,1,subj) ...
        health(:,2,subj) - health(:,1,subj)];
    for t = 1:101 % one regression per time point
        if sum(~isnan(mouseVar(t,:))) > minTrials
            stats = regstats(mouseVar(t,:), rateDiff, 'linear', {'tstat'});
            th.coef(t,1:2,subj) = stats.tstat.beta(2:3);
            th.pval(t,1:2,subj) = stats.tstat.pval(2:3);
        end
    end
    % st calc: get one-sided t-test
    test = th.pval(:,1:2,subj)./2 < .05;
    th.st(subj,1) = 1+find(test(:,1)==0,1,'last');
    th.st(subj,2) = 1+find(test(:,2)==0,1,'last');
    
end
% omit significance times of t >= post-choice
% if sig time >= 101, then it's not really a signif. time, is it?
l.st(l.st>=101) = NaN;
th.st(th.st>=101) = NaN;

% now group significance times
% all subjects
l.stAll = 1+find(ttest(l.coef',0,.05,'right')==0,1,'last');
th.stAll(1) = 1+find(ttest(squeeze(th.coef(:,1,:))',0,.05,'right')==0,1,'last');
th.stAll(2) = 1+find(ttest(squeeze(th.coef(:,2,:))',0,.05,'right')==0,1,'last');
% better SC
ind = successratio > nanmean(successratio);
th.stBetter(1) = 1+find(ttest(squeeze(th.coef(:,1,ind))',0,.05,'right')==0,1,'last');
th.stBetter(2) = 1+find(ttest(squeeze(th.coef(:,2,ind))',0,.05,'right')==0,1,'last');
% worse SC
ind=successratio < nanmean(successratio);
th.stWorse(1) = 1+find(ttest(squeeze(th.coef(:,1,ind))',0,.05,'right')==0,1,'last');
th.stWorse(2) = 1+find(ttest(squeeze(th.coef(:,2,ind))',0,.05,'right')==0,1,'last');


%% time to reach p% full height

% coef. to use:
coef = th.coef;
nRegs = size(coef,2);

percents = .01:.01:.9;
percentsSq = percents.^2;
percentsCub = percents.^3;

timeP = NaN(length(percents),2,28);
intercept = NaN(2,28);
yPredict = NaN(length(percents)+2,2,28);

for subj = 1:28
    for nR = 1:nRegs

        if nRegs> 1
            theseCoef = coef(:,nR,subj);
        else
            theseCoef = coef(:,subj);
        end
        
        % find time at which reach a certain percent
        for p = 1:length(percents)
            target = percents(p) * theseCoef(end); % x percent of final coef.
            
            aboveTarget = (target > 0 & theseCoef >= target) | ... % what times have reached this
                (target < 0 & theseCoef <= target); % call above target if final coef. negative... why?
            
            % if all in this time period are above target:
            timeP(p,nR,subj) = find(~aboveTarget,1,'last')+1;
            
        end
        timeP(timeP==102) = NaN; % never above target.
        
        % fit model
        if any(~isnan(timeP(:,nR,subj)))            
            stats = regstats(timeP(:,nR,subj), ...
                [percents' percentsSq' percentsCub']);
            intercept(nR,subj) = stats.beta(1);
            yPredict(:,nR,subj) = stats.beta(1) + ...
                stats.beta(2).*[0 percents 1]' + stats.beta(3).*[0 percentsSq 1]' + ...
                stats.beta(4).*[0 percentsCub 1]';
        end
        
    end
end


%% figure 2B: example trajectories

fNo=1;
subj=1;
for trial = [11, 28]
    
    % stimuli
    left_box = [162 918.5 510 656.5];
    right_box = [1170 918.5 1518 656.5];
    start_box_height = 65;
    start_box_width = 95;
    start_box = [screenDims(1)*.5-start_box_width screenDims(2)*.8-start_box_height ...
        screenDims(1)*.5+start_box_width screenDims(2)*.8+start_box_height];
    start_box([2 4]) = screenDims(2) - start_box([2 4]); %invert

    % mouse trajectory
    xCoordsRaw = xPos(trial,1:nTPs(trial,subj),subj);
    yCoordsRaw = yPos(trial,1:nTPs(trial,subj),subj);
    figure(fNo); hold on; 
    plot([left_box(1) left_box(3)],[left_box(4) left_box(4)],'k')
    plot([left_box(1) left_box(3)],[left_box(2) left_box(2)],'k')
    plot([left_box(1) left_box(1)],[left_box(2) left_box(4)],'k')
    plot([left_box(3) left_box(3)],[left_box(2) left_box(4)],'k')
    plot([right_box(1) right_box(3)],[right_box(4) right_box(4)],'k')
    plot([right_box(1) right_box(3)],[right_box(2) right_box(2)],'k')
    plot([right_box(1) right_box(1)],[right_box(2) right_box(4)],'k')
    plot([right_box(3) right_box(3)],[right_box(2) right_box(4)],'k')

    plot([start_box(1) start_box(3)],[start_box(4) start_box(4)],'k')
    plot([start_box(1) start_box(3)],[start_box(2) start_box(2)],'k')
    plot([start_box(1) start_box(1)],[start_box(2) start_box(4)],'k')
    plot([start_box(3) start_box(3)],[start_box(2) start_box(4)],'k')

    plot(xCoordsRaw,yCoordsRaw,'color','k','linewidth',2);
    xlim([0 screenDims(1)])
    ylim([0 screenDims(2)])
    xlabel('X Coord. (Pixels)','fontsize',20)
    ylabel('Y Coord. (Pixels)','fontsize',20)
    box on
    set(gca,'fontsize',20)
    fNo=fNo+1;
    
end

%% figure 2C: mean paths

subj=1;
trajAllTrialsL = NaN(101,2,280);
trajAllTrialsR = NaN(101,2,280);
for trial = find(~isnan(choice(:,subj)))' % trial
        if choice(trial,subj)==1
             trajAllTrialsL(1:101,1:2,trial) = spatInterp{subj}(:,:,trial);
        elseif choice(trial,subj)==2
             trajAllTrialsR(1:101,1:2,trial) = spatInterp{subj}(:,:,trial);
        end
end

figure(fNo); hold on
title(['Subject ' num2str(subj) ', Mean Paths'],'fontsize',24)
plot([left_box(1) left_box(3)],[left_box(4) left_box(4)],'k')
plot([left_box(1) left_box(3)],[left_box(2) left_box(2)],'k')
plot([left_box(1) left_box(1)],[left_box(2) left_box(4)],'k')
plot([left_box(3) left_box(3)],[left_box(2) left_box(4)],'k')
plot([right_box(1) right_box(3)],[right_box(4) right_box(4)],'k')
plot([right_box(1) right_box(3)],[right_box(2) right_box(2)],'k')
plot([right_box(1) right_box(1)],[right_box(2) right_box(4)],'k')
plot([right_box(3) right_box(3)],[right_box(2) right_box(4)],'k')

plot([start_box(1) start_box(3)],[start_box(4) start_box(4)],'k')
plot([start_box(1) start_box(3)],[start_box(2) start_box(2)],'k')
plot([start_box(1) start_box(1)],[start_box(2) start_box(4)],'k')
plot([start_box(3) start_box(3)],[start_box(2) start_box(4)],'k')

plotTwoDimError(nanmean(trajAllTrialsL(:,1,:),3),nanmean(trajAllTrialsL(:,2,:),3),...
    nanstd(trajAllTrialsL(:,1,:),0,3)./sqrt(sum(choice(:,subj)==1)),...
    nanstd(trajAllTrialsL(:,2,:),0,3)./sqrt(sum(choice(:,subj)==1)),...
    {'lineColor','k','lineWidth',1,'plotType','shaded','patchSaturation',.3,'overlayLines',false});
plotTwoDimError(nanmean(trajAllTrialsR(:,1,:),3),nanmean(trajAllTrialsR(:,2,:),3),...
    nanstd(trajAllTrialsR(:,1,:),0,3)./sqrt(sum(choice(:,subj)==2)),...
    nanstd(trajAllTrialsR(:,2,:),0,3)./sqrt(sum(choice(:,subj)==2)),...
    {'lineColor','k','lineWidth',1,'plotType','shaded','patchSaturation',.3,'overlayLines',false});

box on
xlim([0 screenDims(1)])
ylim([0 screenDims(2)])
xlabel('X Coord. (Pixels)','fontsize',20)
ylabel('Y Coord. (Pixels)','fontsize',20)

box on
xticks = [0 800 1600];
set(gca,'fontsize',20)
fNo = fNo+1;

%% figure 3A: psychometric curve

% psychometric curve
for subj = 1:28
    
    % mouse tracking trials
    leftChoice = choice(:,subj)==1;
    diffL = like(:,1,subj) - like(:,2,subj);
    total_foods = [sum(diffL==-4); sum(diffL==-3); sum(diffL==-2); ...
        sum(diffL==-1); sum(diffL==0); sum(diffL==1); sum(diffL==2); ...
        sum(diffL==3); sum(diffL==4)];
    left_choice = [sum(diffL==-4 & leftChoice); ...
        sum(diffL==-3 & leftChoice); ...
        sum(diffL==-2 & leftChoice); ...
        sum(diffL==-1 & leftChoice); ...
        sum(diffL==0 & leftChoice); ...
        sum(diffL==1 & leftChoice); ...
        sum(diffL==2 & leftChoice); ...
        sum(diffL==3 & leftChoice); ...
        sum(diffL==4 & leftChoice)];
    prob_left_m(:,subj) =  left_choice ./ total_foods;
    % individual logits:
    b = mnrfit(diffL,[ones(length(leftChoice),1), leftChoice]);
    likeBeta.mouse(subj) = b(2);
    
    % keyboard trials
    leftChoice = key_choice(:,subj)==1;
    diffL = key_like(:,1,subj) - key_like(:,2,subj);
    total_foods = [sum(diffL==-4); sum(diffL==-3); sum(diffL==-2); ...
        sum(diffL==-1); sum(diffL==0); sum(diffL==1); sum(diffL==2); ...
        sum(diffL==3); sum(diffL==4)];
    left_choice = [sum(diffL==-4 & leftChoice); ...
        sum(diffL==-3 & leftChoice); ...
        sum(diffL==-2 & leftChoice); ...
        sum(diffL==-1 & leftChoice); ...
        sum(diffL==0 & leftChoice); ...
        sum(diffL==1 & leftChoice); ...
        sum(diffL==2 & leftChoice); ...
        sum(diffL==3 & leftChoice); ...
        sum(diffL==4 & leftChoice)];
    prob_left_k(:,subj) = left_choice ./ total_foods;    
    % individual logits:
    b = mnrfit(diffL,[ones(length(leftChoice),1), leftChoice]);
    likeBeta.key(subj) = b(2);

end
predictorArray = [];
outcomeArray = [];outcomeArray_m = [];outcomeArray_k = [];
for subj=1:28
    predictorArray = [predictorArray; (-4:4)'];
    outcomeArray_m = [outcomeArray_m; prob_left_m(:,subj)];
    outcomeArray_k = [outcomeArray_k; prob_left_k(:,subj)];
end
% non-linear regression fit:
modelfun = @(b,x)(1 ./ (1+exp(b(1) * x)));
beta0 = 1;
b_left_m = nlinfit(predictorArray,outcomeArray_m,modelfun,beta0);
b_left_k = nlinfit(predictorArray,outcomeArray_k,modelfun,beta0);

figure(fNo); hold on;
colors = [0 0 0; 0 0 0];
% mouse
[~,R,J,CovB,MSE] = nlinfit(predictorArray,outcomeArray_m,modelfun,beta0);
[yPredict,delta] = nlpredci(modelfun,-4:4,b_left_m,R,'Covar',CovB);

temp=shadedErrorBar(-4:4,100*modelfun(b_left_m,-4:4),100*delta,...
    {'color',colors(1,:),'linewidth',1.5},1);
key(1)=temp.mainLine;
% key
[~,R,J,CovB,MSE] = nlinfit(predictorArray,outcomeArray_k,modelfun,beta0);
[yPredict,delta] = nlpredci(modelfun,-4:4,b_left_k,R,'Covar',CovB);
temp=shadedErrorBar(-4:4,100*modelfun(b_left_k,-4:4),100*delta,...
    {'color',colors(2,:),'linewidth',1.5,'linestyle','--'},1);
key(2)=temp.mainLine;
set(gca,'XTick',-4:4)
set(gca,'YTick',0:20:100)
xlim([-4 4])
ylim([0 100])
legend(key,{'mouse','keyboard'},'location','northwest','fontsize',20)
legend BOXOFF
ylabel('% Choose Left','fontsize',24);
xlabel('Left Liking - Right Liking','fontsize',24)
set(gca,'fontsize',20)
fNo = fNo+1;


%% figure 3B: RT by difficulty
% subj 21 has several outlier RTs (was asking a question of experimenter)

mouse_ABS=NaN(5,28);
keyboard_ABS=NaN(5,28);
for subj = 1:28
    
    diff = abs(like(:,1,subj) - like(:,2,subj));
    key_diff = abs(key_like(:,1,subj) - key_like(:,2,subj));
    
    diffArray = 0:4;
    for d=1:length(diffArray)
        mouse_ABS(d,subj) = nanmean(RT(diff==diffArray(d) & ~exclude(:,subj), subj));
        if subj ~=21 || (subj == 21 && diffArray(d) ~= 4)% to exclude crazy outlier
            keyboard_ABS(d,subj) = nanmean(key_RT(key_diff==diffArray(d),subj));
        end
    end
    
    % mouse: regressions
    temp = regstats(RT(~exclude(:,subj),subj), diff(~exclude(:,subj)), 'linear');
    RT_mouse_linreg(subj) = temp.beta(2);
    
    % key:
    temp = regstats(key_RT(key_RT(:,subj)<=2.5,subj), ...
        key_diff(key_RT(:,subj)<=2.5), 'linear');
    RT_key_linreg(subj) = temp.beta(2);
    
end
 
colors = [0 0 0; 0 0 0];
clear key
figure(fNo);hold on;
temp=shadedErrorBar(0:4,nanmean(mouse_ABS.*1000,2),nanstd(mouse_ABS.*1000,0,2)./sqrt(28),...
    {'color',colors(1,:),'linewidth',1.5},1);
key(1)=temp.mainLine;
temp=shadedErrorBar(0:4,nanmean(keyboard_ABS.*1000,2),nanstd(keyboard_ABS.*1000,0,2)./sqrt(28),...
    {'color',colors(2,:),'linewidth',1.5,'linestyle','--'},1);
key(2)=temp.mainLine;
set(gca,'XTick',0:4)
xlim([-.5 4.5])
xlabel('Difficulty (liking difference)','fontsize',24)
ylabel('Mean (ms)','fontsize',24)
set(gca, 'fontsize', 20)
fNo = fNo+1;

%% figure 3C: liking timecourse regression graph (must run timepoint_regressions first)

lColor=[0 0 0];
timeArray = 1:101;%timeArray = 1:320;
xLimitToPlot = timeArray(end);
figure(fNo); hold on;
line([timeArray(1) timeArray(end)],[0 0],'color',[.4 .4 .4]) % zero-line
shadedErrorBar(timeArray,nanmean(l.coef,2), ...
    nanstd(l.coef')'/sqrt(28), {'Color',lColor,'linewidth',1.5},1);
xlim([timeArray(1) timeArray(end)])
xlabel('Time normalized to 101 points','fontsize',24);
ylabel('Beta','fontsize',24)
set(gca,'XTick',0:20:100,'fontsize',20)
fNo = fNo+1;

%% figure 3 D: control ratio histogram


figure(fNo); hold on;
hist(successratio);
h = get (gca, 'children');
set(h,'FaceColor',[.6 .6 .6],'EdgeColor','k','LineStyle','-')
xlabel('SCSR','fontsize',24);
ylabel('Frequency','fontsize',24);
set(gca,'fontsize',20)
fNo = fNo+1;


%% figure 4A: H & T TP regs

labels = {'taste' 'health'};
linestyles = {'-','-.'};
bwColors={[.3 .3 .3] [.3 .3 .3]};%bwColors={[1 0 0] [0 0 1]};
% plot regression:
n_regs = 2;
figure(fNo); hold on;
line([timeArray(1) timeArray(end)],[0 0],'color',[.4 .4 .4]) % zero-line
for r = 1:2
    temp = shadedErrorBar(timeArray,nanmean(th.coef(:,r,:),3), ...
        nanstd(th.coef(:,r,:),0,3)/sqrt(28), ...
        {'Color',bwColors{r},'linewidth',2,'linestyle',linestyles{r}},1);
    graphkey(r) = temp.mainLine;
end
for r = 1:2
    line([th.stAll(r) th.stAll(r)],...
        get(gca,'YLim'),'color',bwColors{r},'linestyle',linestyles{r},...
        'linewidth',2)
end
xlim([1 101]);ylim([-2 16])
set(gca,'fontsize',20)
legend(graphkey,labels,'location','northwest','fontsize',20);
legend BOXOFF
xlabel('Time normalized to 101 points','fontsize',24);
ylabel('Beta','fontsize',24)
fNo = fNo+1;


%% figure 4B: H T ST dist

mx = nanmax(nanmax(th.st));
mn = nanmin(nanmin(th.st));
xvalues = linspace(mn,mx,10);        

figure(fNo); hold on;
hist(th.st(:,1),xvalues);
hist(th.st(:,2),xvalues);
h = get (gca, 'children');
set(h(2),'FaceColor',[.2 .2 .2],'EdgeColor','k')
set(h(1),'FaceColor',[.6 .6 .6],'EdgeColor','k','LineStyle','-.')
alpha(h(1),.6)
graphkey = h;
legend(graphkey,{'health','taste'},'location','northwest','fontsize',20)
legend BOXOFF
xlim([40,105])
xlabel('Significance time','fontsize',24);
ylabel('Frequency','fontsize',24);
set(gca,'fontsize',20)
fNo = fNo+1;


%% figure 5 A: st by SC

% time course regression
clear labels
catCols{1} = {[.8 0 0],[.8 0 0]};% for taste
catCols{2} = {[0 0 .8],[0 0 .8]};% for health
linestyles{1} = {'-','--'}; % for taste
linestyles{2} = {'-','--'}; % for health
labels{1} = 'taste, better control';
labels{2} = 'taste, worse control';
labels{3} = 'health, better control';
labels{4} = 'health, worse control';
subjTypeNames = {'stBetter','stWorse'};
clear graphkey; graphkey_count = 1;
clear coefToPlotBetter coefToPlotWorse
figure(fNo); hold on;
line([timeArray(1) timeArray(end)],[0 0],'color',[.6 .6 .6]) % zero-line
for r = 1:2
    for subjType = 1:2
        coefToPlot = nanmean(th.coef(:,r,controlInd{subjType}),3);
        stdForPlot = (nanstd(th.coef(:,r,controlInd{subjType}),0,3)/sqrt(sum(controlInd{subjType})))';
        temp = shadedErrorBar(timeArray,coefToPlot, stdForPlot, ...
            {'Color',catCols{r}{subjType},'linewidth',2.5,'linestyle',linestyles{r}{subjType}},1);
        graphkey(graphkey_count) = temp.mainLine;
        graphkey_count = graphkey_count + 1;
    end
end
yLimits = get(gca,'YLim');
for r = 1:n_regs
    for subjType = 1:2
        line([th.(subjTypeNames{subjType})(r) ...
            th.(subjTypeNames{subjType})(r)],...
            yLimits,'color',catCols{r}{subjType},...
            'linestyle',linestyles{r}{subjType},'linewidth',2.5);
    end
end
xlim([timeArray(1) timeArray(end)+1])
ylim([-2 17])
set(gca,'fontsize',20)
xlabel('Time normalized to 101 points','fontsize',24);
ylabel('Beta','fontsize',24)
fNo = fNo+1;


%% figure 5 B

% significance time diff by control success
colorToUse='k';
clear graphkey
figure(fNo); hold on;
plot(th.st(:,1)-th.st(:,2),...
     successratio,'.','color',colorToUse,'markersize',40)
xPlotLocation = get(gca,'XLim');
stats = regstats(successratio,th.st(:,1)-th.st(:,2));
xvals = linspace(xPlotLocation(1),xPlotLocation(2));
yvals = (stats.beta(1)+stats.beta(2).*xvals);
plot(xvals,yvals,'color',colorToUse,'linewidth',1.5);
if sum(isnan(th.st(:,1)-th.st(:,2))) > 0
    graphkey = plot(ones(sum(isnan(th.st(:,1)-th.st(:,2))),1)*xPlotLocation(1),...
        successratio(isnan(th.st(:,1)-th.st(:,2))),...
        'marker','x','color',colorToUse,'markersize',15,'linewidth',2,'linestyle','none');
end
xlim([xvals(1) xvals(end)])
ylabel('SCSR','fontsize',24)
xlabel({'Significance time: Taste - health'},'fontsize',24)
set(gca,'fontsize',20)
fNo = fNo+1;


%% figure 6A: percent max

figure(fNo); hold on;
linestyles = {'-','-.'};
clear graphkey
bwColors={[.3 .3 .3] [.3 .3 .3]};%bwColors={[1 0 0] [0 0 1]};
for r = 1:2
    y = nanmean(timeP(:,r,:),3);
    stdForPlot = nanstd(timeP(:,r,:),0,3)'/sqrt(28)';
    temp = shadedErrorBar(percents,y,stdForPlot, ...
        {'Color',bwColors{r},'linewidth',2.5,'linestyle',linestyles{r}},1);
    graphkey(r) = temp.mainLine;
end
legend(graphkey,{'taste','health'},'location','northwest','fontsize',20)
legend BOXOFF
xlabel('Proportion of full attribute weighting','fontsize',24)
ylabel('Normalized time','fontsize',24)
set(gca,'fontsize',20)
xlim([0 1]); ylim([35 95])
fNo = fNo+1;


%% figure 6B

y=th.st(:,1)-th.st(:,2);
x=intercept(1,:) - intercept(2,:);

figure(fNo); hold on;
plot(x,y,'ko','markersize',15,'linewidth',2)
xPlotLocation = get(gca,'XLim');
plot(x(isnan(y)),ones(sum(isnan(y)),1)*-50,'xk','markersize',15,'linewidth',2)
stats = regstats(y,x);
xvals = linspace(xPlotLocation(1),xPlotLocation(2));
yvals = (stats.beta(1)+stats.beta(2).*xvals);
plot(xvals,yvals,'color','k','linewidth',1.5);
xlim(xPlotLocation)
ylim([-50 40])
ylabel({'Significance time:','Taste-health'},'fontsize',24)
xlabel({'Intercept: Taste-health'},'fontsize',24)
set(gca,'fontsize',20)
fNo = fNo+1;


%% figure 6C

y=successratio;
x=intercept(1,:) - intercept(2,:);

figure(fNo); hold on;
plot(x,y,'k^','markersize',15,'linewidth',2,'markerfacecolor','k')
stats = regstats(y,x);
xvals = linspace(xPlotLocation(1),xPlotLocation(2));
yvals = (stats.beta(1)+stats.beta(2).*xvals);
plot(xvals,yvals,'color','k','linewidth',1.5);
xlim([xvals(1) xvals(end)])
ylabel({'SCSR'},'fontsize',24)
xlabel({'Intercept: Taste-health'},'fontsize',24)
set(gca,'fontsize',20)
fNo = fNo+1;


%% figure 7A: attribute weighting and SCSR

for subj=1:28
     tD = taste(:,2,subj)-taste(:,1,subj);
     hD = health(:,2,subj)-health(:,1,subj);
     x = [tD hD];
     y = choice(:,subj) == 2;
     [~,~,stats]=glmfit(x,y,'binomial','logit');
     tW(subj) = stats.beta(2);
     hW(subj) = stats.beta(3);
end
tW(tW>10) = NaN;
% plot it
thisHcolor=[0 0 0]; %thisHcolor=[0 0 1];
thisTcolor=[0 0 0]; %thisTcolor=[1 0 0];
figure(fNo); hold on;
plot(hW,successratio,'marker','s','linestyle','none','color',thisHcolor, ...
    'markerfacecolor',thisHcolor, 'markersize',15)
plot(tW,successratio,'marker','o','linestyle','none','color',thisTcolor, ...
    'linewidth',2, 'markersize',15)
xPlotLocation = get(gca,'XLim');
stats = regstats(successratio,tW);
xvals = linspace(xPlotLocation(1),xPlotLocation(2));
yvals = (stats.beta(1)+stats.beta(2).*xvals);
plot(xvals,yvals,'color',thisTcolor,'linewidth',1.5);
stats = regstats(successratio,hW);
xvals = linspace(xPlotLocation(1),xPlotLocation(2));
yvals = (stats.beta(1)+stats.beta(2).*xvals);
plot(xvals,yvals,'color',thisHcolor,'linewidth',1.5,'linestyle','--');
key(1) = plot(-500,-500,'marker','s','linestyle','--','color',thisHcolor,'linewidth',2,'markerfacecolor',thisHcolor,'markersize',10);
key(2) = plot(-500,-500,'marker','o','linestyle','-','color',thisTcolor,'linewidth',2,'markersize',10);
xlim(xPlotLocation)
ylim([0 max(successratio)])
set(gca,'fontsize',20)
legend(key,{'Health','Taste'},'fontsize',20)
legend BOXOFF
xlabel('Attribute weighting on choice','fontsize',24)
ylabel('SCSR','fontsize',24)
fNo = fNo+1;


%% fig 7B: decision weights by ST

y=tW-hW;
x=th.st(:,1)-th.st(:,2);
figure(fNo); hold on;
plot(x,y,'k.','markersize',40)
plot(ones(sum(isnan(x)),1)*-50,y(isnan(x)),'xk','markersize',15,'linewidth',2)
xPlotLocation = get(gca,'XLim');
set(gca,'fontsize',20)
stats = regstats(y,x);
xvals = linspace(xPlotLocation(1),xPlotLocation(2));
yvals = (stats.beta(1)+stats.beta(2).*xvals);
plot(xvals,yvals,'color','k','linewidth',1.5);
xlim(xPlotLocation)
xlabel({'Significance time: Taste-health'},'fontsize',24)
ylabel({'Attribute weighting:','Taste-health'},'fontsize',24)
fNo = fNo+1;



%% stats in abstract

% normalized time diff b/n taste and health ST, multiplied by RT
for subj = 1:length(subjects)
    % get mean in milliseconds
     mtRT(subj) = nanmean(RT(:,subj)*1000);
end
aveST = nanmean(th.st)/101; % percent of trial each ST represents
ST_perc_Diff = (aveST(2) - aveST(1)); % percent of trial difference
tasteEarlierMS = nanmean(mtRT) * ST_perc_Diff;
sprintf(['"taste attributes are processed about %0.0f ms earlier than '...
    'health attributes during the choice process"'],...
    tasteEarlierMS)

% percent variance explained
% omit subjects with no ST for health:
stats = regstats(successratio, th.st(:,1)-th.st(:,2));
percOne = stats.rsquare * 100;
% for subjects without an ST, replace that value with 102 (one timepoint
% past the final timepoint):
t = th.st(:,1);
t(isnan(t))= 102;
h = th.st(:,2);
h(isnan(h)) = 102;
stats = regstats(successratio,t-h);
percTwo = stats.rsquare * 100;
sprintf(['"We also find that %0.0f - %0.0f%% of observed individual differences '...
    'in self-control ability can be explained by differences in the relative '...
    'speed with which taste and health attributes are processed."'], ...
    percOne,percTwo)

%% methods: cursor x velocity

for subj = 1:length(subjects)
    % velocity:
    for nT = find(~exclude(1:nTrials(subj),subj))'
        veloc(:,nT) = temporal{subj}(2:end,1,nT) - temporal{subj}(1:end-1,1,nT);
    end
    % test: velocity > 0?
    for t=1:100
        [~,p(t,subj)] = ttest(abs(veloc(t,:)),0,'tail','right');
    end
end
time = find(sum(p<.01,2) == 28,1,'last');
sprintf(['Cursor velocity was significantly greater than zero for normalized '...
    'time points 1 through %d (two-sided t-test; p<0.01, uncorrected for '...
    'multiple comparisons), suggesting that movements were fluid on average.'],...
    time)

%% methods: trial exclusion

for subj = 1:length(subjects)
    percExc(subj) = sum(exclude(1:nTrials(subj),subj)) ./ nTrials(subj);
    percYCross(subj) = sum(reason.yCross(1:nTrials(subj),subj)) ./ nTrials(subj);
end

% trial exclusions:
sprintf(['...excluded from further analysis (mean %0.1f%% of trials).'], ...
    nanmean(percExc*100))

sprintf(['For comparison, the mean RT in mouse trials was %0.0f ms (S.D. = '...
    '%0.0f ms).'], nanmean(nanmean(RT*1000)),nanstd(nanmean(RT*1000)))

sprintf(['Second, we also removed a small fraction of trials in which the '...
    'mouse trajectory crossed the y-axis more than three times (mean %0.1f%% '...
    'of trials per subject; %d of 28 subjects had at least one such trial).'],...
    nanmean(percYCross*100),sum(percYCross>0))
    

%% results: paradigm validation

% choice curve
[~,p,~,stats]=ttest(likeBeta.mouse,likeBeta.key);
sprintf(['We found that the two cases led to indistinguishable choice '...
    'curves (mean logistic slope mouse = %0.2f, mean logistic slope '...
    'keyboard = %0.2f; p=%0.2f, t(%0.0f)=%0.2f, two-tailed).'], ...
    nanmean(likeBeta.mouse),nanmean(likeBeta.key),p,stats.df,stats.tstat)

% reaction times
[~,p,~,stats]=ttest(nanmean(RT),nanmean(key_RT));
sprintf(['Although mean reaction times were longer in the mouse condition '...
    '(mouse mean = %0.0f ms; keyboard mean = %0.0f ms; t(%0.0f)=%0.2f, p=%0.4f, '...
    'two-tailed)'],nanmean(nanmean(RT*1000)), nanmean(nanmean(key_RT*1000)),...
    stats.df, stats.tstat, p)

% reaction time and choice difficulty
[~,p,~,stats]=ttest(RT_mouse_linreg,RT_key_linreg);
sprintf(['the absolute value of the liking rating differences (mouse mean '...
    'linear regression slope = %0.03f; keyboard mean linear regression '...
    'slope = %0.3f; p=%0.2f, t(%0.0f)=%0.2f, two-tailed).'],...
    nanmean(RT_mouse_linreg), nanmean(RT_key_linreg), p, stats.df, stats.tstat)
    
% time to significance for liking value rating
sprintf(['We found that, across the group, the earliest time at which the '...
    'value information had a sustained, significant, effect on the '...
    'trajectories, without ever returning to being not-significant, for the '...
    'group was t = %0.0f.'],l.stAll)

% dietary self-control success ratio (SCSR)
sprintf(['SCSRs across the group, and shows that there was substantial '...
    'variation across subjects (mean = %0.0f%%, range: %0.0f-%0.0f%%).'], ...
    nanmean(successratio*100), min(successratio*100), max(successratio*100))

% reaction times in control-requited and not-required cases
[~,p,~,stats]=ttest(rtSCReq, rtSCNotReq);
sprintf(['There were no significant RT differences between challenge and '...
    'non-challenge trials (challenge mean = %0.0f ms; non-challenge mean = '...
    '%0.0f ms; t(%0.0f)=%0.2f, p=%0.2f, two-tailed)'], nanmean(rtSCReq*1000), ...
    nanmean(rtSCNotReq*1000), stats.df, stats.tstat, p)


%% results: Taste is reflected in the choice process earlier than health

sprintf(['mouse trajectories were influenced by taste at significantly '...
    'earlier time slices compared to health (taste: t = %0.0f; health: '...
    't = %0.0f; p<0.05 threshold, one-tailed).'], th.stAll)

noTasteHealthST = (sum(isnan(th.st)));
sprintf(['We found that health never became a significant influence on the '...
    'trajectory for %0.0f%% of subjects, whereas for taste this happened for '...
    '%0.0f subjects'],100*(noTasteHealthST(2)/28), noTasteHealthST(1))

[~,p,~,stats]=ttest(th.st(:,1), th.st(:,2),'tail','left');
sprintf(['the mean earliest time for a significant effect of health was '...
    'later than the mean earliest time for taste (mean difference = %0.2f, '...
    'p=%0.04f, t(%0.0f)=%0.2f, one-tailed).'], nanmean(th.st(:,2) - th.st(:,1)), ...
    p, stats.df, stats.tstat)

sprintf(['taste information began to influence the choice process about '...
    '%0.0f%% earlier than health information.'], ...
    ((th.stAll(2) - th.stAll(1)) / 101)*100)

sprintf(['onset of health attribute processing by the choice circuitry '...
    'approximately %0.0f ms later than the onset for taste.'],tasteEarlierMS)


%% results: Individual differences in dietary self-control are associated with differences in the relative speed at which taste and health attributes are computed.

sprintf(['in the high self-control group, the paths for taste and health '...
    'were quite similar and became significantly greater than zero at '...
    'approximately the same time (health: t = %0.0f; taste: t = %0.0f)'], ...
    fliplr(th.stBetter))

sprintf(['for the low self-control group, the latency at which taste became '...
    'significant was similar to that of the high self-control group '...
    '(t = %0.0f), but that for health occurred much later (t = %0.0f).'], ...
    th.stWorse)

[~,p,~,stats]=ttest(th.st(controlInd{1},1), th.st(controlInd{1},2),'tail','left');
sprintf(['The normalized time at which taste and health became significant '...
    'did not differ for high self-control subjects (mean health: t = %0.2f; '...
    'mean taste: t = %0.2f; mean difference = %0.2f; t(%0.0f)=%0.2f, p=%0.2f, '...
    'one-tailed)'],fliplr(nanmean(th.st(controlInd{1},:))), ...
    nanmean(th.st(controlInd{1},1) - th.st(controlInd{1},2)), ...
    stats.df, stats.tstat, p)

[~,p,~,stats]=ttest(th.st(controlInd{2},1), th.st(controlInd{2},2),'tail','left');
sprintf(['but they were significantly different in the low '...
    'self-control group (mean health: t = %0.2f; mean taste: t = %0.2f; '...
    'mean difference = %0.2f, t(%0.0f)=%0.2f, p=%0.2f, one-tailed).'], ...
    fliplr(nanmean(th.st(controlInd{2},:))), ...
    nanmean(th.st(controlInd{2},1) - th.st(controlInd{2},2)), ...
    stats.df, stats.tstat, p)
    
% time to attribute significance vs. SCSR
t=th.st(:,1);
t(isnan(t))= 102;
h=th.st(:,2);
h(isnan(h)) = 102;
stats = regstats(successratio,t-h);
sprintf(['In one case we assumed that, for the problem subjects, health '...
    'first became significant at t = 102, one time point after the final '...
    'time point, and estimated the regression using all of the subjects '...
    '(linear regression slope = %0.03f, p=%0.04f, R^2=%0.02f).'], ...
    stats.beta(2), stats.fstat.pval, stats.rsquare)

percNoTasteHealthST = sum(isnan(th.st)) ./ 28;
stats = regstats(successratio,th.st(:,1)-th.st(:,2));
sprintf(['In the other case, we simply excluded these subjects, carrying out '...
    'the regression using only the %0.00f%% of subjects for whom health became '...
    'significant during the trial. This yielded a similar result (linear '...
    'regression slope = %0.03f, p=%0.04f, R^2=%0.02f).'], ...
    (1-percNoTasteHealthST(2))*100, stats.beta(2), stats.fstat.pval, stats.rsquare)

% time to reach p% height: correlation with signif. time:
stats = regstats(th.st(:,1)-th.st(:,2), intercept(1,:) - intercept(2,:));
sprintf(['the differences in earliest significant time identified with the '...
    'previous analysis are highly correlated with those identified with this '...
    'new method (R^2=%0.2f, p=%0.5f).'], stats.rsquare, stats.fstat.pval)

% time to reach p% height: correlation with SCSSR
stats = regstats(successratio,intercept(1,:) - intercept(2,:));
sprintf(['Likewise, the alternative measure of the relative earliest times '...
    'can still significantly predict a sizable fraction of the individual '...
    'differences in self-control (Fig. 6C; R^2=%0.2f, p=%0.5f).'], ...
    stats.rsquare, stats.fstat.pval)

sprintf(['Together, these results suggest that a sizable fraction of '...
    'individual differences in dietary self-control, on the order of '...
    '%0.0f - %0.0f%%'], percOne,percTwo)


%% results: Relationship between relative computation time and attribute decision weights

% weight of health and taste on choice
stats = regstats(successratio,tW);
tasteP = stats.fstat.pval;
tasteRsq = stats.rsquare;
stats = regstats(successratio,hW);
healhP = stats.fstat.pval;
healthRsq = stats.rsquare;
sprintf(['Both attributes were highly significant and, on their own, '...
    'explained a sizable fraction of the individual differences in '...
    'self-control (health: p=%0.4f, R^2=%0.2f; taste: p=%0.4f, R^2=%0.4f).'], ...
    healhP, healthRsq, tasteP, tasteRsq)

% attribute weighting in choice vs. significance time
stats = regstats(tW-hW, th.st(:,1)-th.st(:,2));
sprintf(['As shown in Figure 7B, we found that the two variables had a '...
    'significant and sizable relation (slope = %0.2f, p=%0.3f, R^2=%0.4f).'], ...
    stats.beta(2), stats.fstat.pval, stats.rsquare)

% mediation analysis
% taste
nR=1;
IV = th.st(:,nR);
DV = successratio';
med = squeeze(th.coef(end, nR, :));
% 1. IV predicts DV?
stats = regstats(DV,IV);
stepOne=[stats.beta(2) stats.tstat.pval(2)];
% 2. mediator predicted by IV?
stats = regstats(med,IV);
stepTwo=[stats.beta(2) stats.tstat.pval(2)];
    % 2b. mediator predicts DV? (incidental for med. analysis, but use in figure)
    stats = regstats(DV,med);
    stepTwoB=[stats.beta(2) stats.tstat.pval(2)];
% 3. both in predictor model
stats = regstats(DV,[IV med]);
stepThree = [stats.beta(2:3) stats.tstat.pval(2:3)]; %[betas column, p column]
% 4. percent change in path strength:
changePathStrngth=(abs(stepOne(1)) - abs(stepThree(1,1))) / abs(stepOne(1)); 
sprintf(['Using this method, we found that the final weight of taste was a '...
    'significant mediator of taste’s time-to-significance, reducing path '...
    'strength by %0.2f%%, and making taste’s time-to-significance '...
    'not-significant for predicting SCSR (p=%0.2f).'], changePathStrngth*100,...
    stepThree(1,2))

% health
nR=2;
IV = th.st(:,nR);
DV = successratio';
med = squeeze(th.coef(end, nR, :));
% 1. IV predicts DV?
stats = regstats(DV,IV);
stepOne=[stats.beta(2) stats.tstat.pval(2)];
% 2. mediator predicted by IV?
stats = regstats(med,IV);
stepTwo=[stats.beta(2) stats.tstat.pval(2)];
    % 2b. mediator predicts DV? (incidental for med. analysis, but use in figure)
    stats = regstats(DV,med);
    stepTwoB=[stats.beta(2) stats.tstat.pval(2)];
% 3. both in predictor model
stats = regstats(DV,[IV med]);
stepThree = [stats.beta(2:3) stats.tstat.pval(2:3)]; %[betas column, p column]
% 4. percent change in path strength:
changePathStrngth=(abs(stepOne(1)) - abs(stepThree(1,1))) / abs(stepOne(1)); 
sprintf(['This was also true for the health coefficient, reducing path '...
    'strength by %0.2f%%, and also making health’s time-to-significance '...
    'non-significant for predicting SCSR (p=%0.2f).'], changePathStrngth*100,...
    stepThree(1,2))



end