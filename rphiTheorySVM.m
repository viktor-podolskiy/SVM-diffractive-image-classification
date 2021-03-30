% SVM training and evaluation script 
%
% part of the diffractive image classification with subwavelength resolution project
% (c) 2021 V. Podolskiy - University of Massachusetts Lowell

% see A. Ghosh, et.al., ACS Photonics (2021) for more details 

clear 
rng(1111)

ln=5; % the value of the "l" parameter from the manuscript

numTrain=10; %number of images to train SVM 
numTest=15; %number of test images to evaluate SVM performance; numTrain+numTest must be smaller than total number of images
numTry=10; %number of independent SVM tries for each {m,j,ln} combination

%setup m,j array
m1Lst=(0:2:20); 
j1Lst=(1:2:20); 
[m2,j2]=meshgrid(m1Lst,j1Lst); 

accuracy=0*m2; %accuracy array, based on SVM diagnostic data
accuracyMan=0*m2; %accuracy array, based on actual testing data
accuracyCat=zeros(numel(accuracyMan),12); 

for it=1:numel(m2)
    m1st=m2(it); 
    j1st=j2(it); 

    m1=(1:ln)-1+m1st; 
    j1=(1:ln)-1+j1st;  
    [fullTbl,labels]=rphi2tableFun(m1,j1,'rphiImageData.mat'); 
    predictLbls=cell(length(labels),1); 
    for ipl=1:length(predictLbls)
        predictLbls{ipl}=labels{ipl}; 
    end 

    if ~exist('categories','var') 
        categories=unique(fullTbl.fullLabels); 
    end 
    for itry=1:numTry
        % split the table into the training table and testing table
        trainTbl=fullTbl([],:); 
        testTbl=fullTbl([],:);
        for icat=1:length(categories)
            subsetTbl=fullTbl(strcmp(fullTbl.fullLabels,categories{icat}),:); 

            indPerm=randperm(size(subsetTbl,1)); 
            trainTbl=[trainTbl;subsetTbl(indPerm(1:numTrain),:)]; 
            testTbl=[testTbl; subsetTbl(indPerm(numTrain+(1:numTest)),:)]; 
        end 


        [classifier,valid]=trainLinearSVM(trainTbl,predictLbls,categories); 
        accuracy(it)=accuracy(it)+valid; 

        % detect
        detectTbl=testTbl(:, (1:size(testTbl,2)-1)); 
        imgPrediction=classifier.predictFcn(detectTbl);
        imgPrediction=table(imgPrediction,'VariableNames',{'detectLabels'});
        tmp=strcmp(testTbl.fullLabels,imgPrediction.detectLabels);
        accuracyMan(it)=accuracyMan(it)+sum(tmp)/length(tmp); 
        for icat=1:length(categories)
            tmp=imgPrediction.detectLabels(strcmp(testTbl.fullLabels,categories{icat})); 
            tmp=strcmp(tmp,categories{icat}); 
            accuracyCat(it,icat)=accuracyCat(it,icat)+sum(tmp)/length(tmp); 
        end 
        
        
    end 
    accuracy(it)=accuracy(it)/numTry; 
    accuracyMan(it)=accuracyMan(it)/numTry; 
    accuracyCat(it,:)=accuracyCat(it,:)/numTry; 

    disp(['iter:',num2str(it),' of ',num2str(numel(m2))]); 
end 

%% save results
save('SVMperformanceSummary.mat',...
    'accuracy','accuracyMan', 'accuracyCat', ...
    'numTrain','numTry','numTest',...
    'm2','j2','ln','categories')

%% plot the summary

% plot overall SVM accuracy - surface plots
figure(10)
clf
subplot(1,2,1)
surf(m2,j2,accuracy*100, 'LineStyle','none')
set(gca,'FontSize',18)
xlabel('m')
ylabel('j')
ylim([1 max(j2(:))])
view(2)
colormap('jet')
set(gca,'CLim',[70. 92])
colorbar

subplot(1,2,2)
surf(m2,j2,accuracyMan*100, 'LineStyle','none')
set(gca,'FontSize',18)
xlabel('m')
ylabel('j')
ylim([1 max(j2(:))])
view(2)
colormap('jet')
set(gca,'CLim',[70. 92])
colorbar


% plot object-specific accuracy - bar plots
indOrder=[6,7,2,5,3,10,8,4,1,12,9,11]';

xlbs=categories; 

figure(4)
clf
indCorrBar=100*(sum(accuracyCat,1)/size(accuracyCat,1)).'; 
bar(indCorrBar(indOrder),'stacked')
set(gca,'FontSize',18)
xtickangle(90)
xticks((1:length(indOrder)))
xticklabels(xlbs(indOrder))
ylabel('Accuracy')


