%% load data and variables
clear all
close all
addpath 'U:\eng_research_nialab\users\suezhang\Dispersion Indices and Autophagy Processing'
%%%%%%%%%%%%%%%%%%  find image folders %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%If I do more than 1 treatment at a time it doesn't separate the data so I
%need to fix that
%analyzeimageStacks -- does this bit depth need to be 16 or 12?
bit_depth=input('What is the bit depth?'); %usually 16 for widefield and 12 for confocal

disp('Select the C1 files')
[filename, pathname] = uigetfile('*.tif','Multiselect','on'); %select the C1 files and it will recognize the other channels
addpath(pathname)
if iscell(filename)
    fname = filename;
elseif filename ~= 0
    fname{1} = filename;   
else
    disp('No files selected')
    return
end

for fileindex=1:length(fname)
saveName = erase(fname{fileindex},{'C1-MAX_','.tif'}); 

pNames = strings(10,length(pathname)) ; 
pNames =saveName;

fileFolder1{1,1} = pathname;

maskPath=pathname;
parent_name=saveName

disp("Locating Images: ")

    
 
       Img(1,1,1) = strcat({pathname},insertBefore(saveName,'SZ','C1-MAX_'),'.tif'); %Nucleus
       Img(2,1,1) = strcat({pathname},insertBefore(saveName,'SZ','C2-MAX_'),'.tif'); %Actin
       Img(3,1,1) = strcat({pathname},insertBefore(saveName,'SZ','C3-MAX_'),'.tif'); %Mitochondria
       Img(4,1,1) = strcat({maskPath},insertBefore(fname{fileindex},'.tif','_CellMask'));
       Img(4,1,1) = insertAfter(Img(4,1,1),'.tif','f');
       
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% process all the images %%%%%%%%%%%%%%%%%%%%%%%
c1 = datetime('now') ;
count = 1;
cellData = [];
CgroupTreatment = [];
CgroupImage = [];
for ff = 1:length(Img(1,1,:)) % treatment / staining group
    % process image groups
    countColumn = 1;
    for dd =  1:length(Img(1,:,1)) % run through all images in stainng group
        if Img(1,dd,ff) ~= "" 
            % create summary stats about images
            [thresholds(count,:),fData(count,:),cData] = thresholdImages(Img(:,dd,ff),pNames,bit_depth) ;                
            % assign frame grouping info
            FgroupTreatment(count,:) = ff ;
            FgroupImage(count,:) = dd ; 
            
            %assign cell grouping infomation

            cellData = [cellData ; cData] ;
            CgroupTreatment = [CgroupTreatment; cData.*0 + ff];
            CgroupImage = [CgroupImage; cData.*0 + dd];
            
            count  = count + 1;
        
            
        end
    end 
end

thresh_data = [thresholds,fData] ;


c2 = datetime ;
disp("Analysis Complete") ;

Results.data = thresh_data ;
Results.cellData = cellData;
Results.FgroupTreatment = FgroupTreatment ; 
Results.FgroupImage = FgroupImage ;
Results.CgroupTreatment = CgroupTreatment ; 
Results.CgroupImage = CgroupImage ;

Results.processInfo = [strcat("Analysis started at ",string(c1)) ;
                       strcat("Analysis finished at ",string(c2)) ] ;

disp(strcat("Analysis started at ",string(c1))) ;
disp(strcat("Analysis finished at ",string(c2))) ;
clear fileindex
save([saveName,'_analysis']) ; 
disp(strcat("The data has been saved as: ", saveName) );

end
%%
%%
%% Functions
function [groupThresh, frameData, cellData] = thresholdImages(im,name,bit_depth)
% read in images
disp(strcat("Processing Image ",name));

tic
image(:,:,:,1) = zStackReader(im{1,1},bit_depth); %nucleus bit depth 12 originally
image(:,:,:,2) = zStackReader(im{2,1},bit_depth); %actin bit depth 12 originally
image(:,:,:,3) = zStackReader(im{3,1},bit_depth); %mitochondria bit depth 12 originally
image(:,:,:,4) = zStackReader(im{4,1},8); %bit depth is different for mask
 
%data.images = im ; 
%% Threshold image stack using otsu method 
% grouped thresholding of stack
gT(1:2) = multithresh(image(:,:,:,3),2);
gT(3) = multithresh(image(:,:,:,3)) ;
gT(4:5) = multithresh(image(:,:,:,2),2);
gT(6) = multithresh(image(:,:,:,2));


if isfile(im{4,1}) == true
    mask = imread(im{4,1});
    mask = double(mask) ;
    % organize masking information
    cellLocationID = unique(mask(mask>0));
    cellData = zeros([length(cellLocationID) 7]).*-1 ;
    figure()
   imshow(mask)
 
    for ff = 1:length(cellLocationID)
        singleCellID = cellLocationID(ff);
        singleCellMask = mask;
        singleCellMask(singleCellMask~=singleCellID) = 0 ;
        cellData(ff,:) = analyzeimageStacks(image,gT,singleCellMask) ;
        %figure()
        imshow(singleCellMask);
        
    end


else
    disp("No Cell Mask Data Found... Processing Frames Only!")
    cellData = [0 0];
end
% get output values
groupThresh  = gT;
fmask = ones(size(image(:,:,1,1))) ;
frameData = analyzeimageStacks(image,gT,fmask) ;

end


function result = analyzeimageStacks(imData,threshVals,cellmask)
% result  = [  MLD,IMcov,Theil_T,Gini, meanLC3, meanActin, nucleiArea ]
%% detemine input type
if class(imData) == 'double'
    imageStacks_pre = imData ;
    clear imData

    %Does bit depth need to be changed here? Maybe not because it was
    %normalized before using zStackReader
elseif class(imData) == 'string'
    disp(strcat("Reprocessing Image",imData(1)))
    imageStacks_pre(:,:,:,1) = zStackReader(imData{1,1},16); %changed from 16
    imageStacks_pre(:,:,:,2) = zStackReader(imData{2,1},16); %changed from 16
    imageStacks_pre(:,:,:,3) = zStackReader(imData{3,1},16); %changed from 16
    clear imData
else
    disp("Unexpected Data Input")
   disp(strcat("Data Input Class: ",class(imData)))
   return
end
%% add cell mask to stacks
nanmask = find(cellmask) ;
% find cell areas based on stacks
imageStacks = zeros(size(imageStacks_pre));
for pp = 1:length(imageStacks_pre(1,1,1,:))
    for kk = 1:length(imageStacks_pre(1,1,:,1))
        slice = imageStacks_pre(:,:,kk,pp);
        slice(~cellmask) = -1;
        imageStacks(:,:,kk,pp) = slice ;
    end
end


%% locate cell area based on BActin stain
Cell_Area_logic = imageStacks(:,:,:,2)>=threshVals(4) ;
nucleiAll = imageStacks(:,:,:,1) ;
nuclei = nucleiAll >=0 ;
nucleiIM = nucleiAll(nuclei);
nT = multithresh(nucleiIM);
nnT = imageStacks(:,:,:,1) < nT ;

mito_all = imageStacks(:,:,:,3) ; 
imshow(mito_all)
mito_cell_all =mito_all(Cell_Area_logic & mito_all < 0.99 & nnT) ; % must be within cell and not saturated and non-nuclear


%clear Cell_Area_logic 

%% calculate gfp /AF488 measures
IMcov_mito = std(mito_cell_all(:),1,'all') / mean(mito_cell_all(:),'all');
N = length(mito_cell_all);
M = mean(mito_cell_all) ;
Theil_T_mito = 1/N * sum( (mito_cell_all ./M).*log(mito_cell_all ./M) ) ;
% Gini copied from Matlab function, 
Weights = mito_cell_all./sum(mito_cell_all) ; % normalize all data
SortedWeights = sort(Weights) ; % sort from low to high
ProportionPixels = [0 (1:N)./N] ;
ProportionValue = [0 , cumsum(SortedWeights')] ; % basically a cdf?
Gini_mito = 1 - sum(ProportionValue(1:end-1)+ProportionValue(2:end))/N ; % calculate Gini Index
clear Weights SortedWeights ProportionPixels ProportionValue M N 
% calculate mean log deviation
MLD_mito = log(mean(mito_cell_all)) - mean(log(mito_cell_all)) ; 

%% calculate control parameters
m_mito = mean(imageStacks(:,:,:,3),'all');
mActin = mean(imageStacks(:,:,:,2),'all');
nucleiArea = sum(imageStacks(:,:,:,1)>nT,'all');

%% paste output
result = [MLD_mito, IMcov_mito, Theil_T_mito, Gini_mito,m_mito, mActin, nucleiArea] ;
% 1.MLD GE(0)
% 2.COV
% 3.Theil_T GE(1)
% 4.Gini
% 5.meanLC3
% 6. meanActin
% 7. nucleiArea

end

