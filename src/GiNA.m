function [data, filesNamesA, matrixA]=GiNA(folder,resizeIndex,method,minApx,maxApx,realLength,keepOrder,showFig,writeT)
%
% GINA Extracts objects from pictures and measures different parameters.
%
% -------------------------------------------------------------------------------------------
% ** A quick test: data=GINA('/directory/picture.JPG',0.4,[100 3 6 6],700,1000,2,false,true,true)
% -------------------------------------------------------------------------------------------
%            data=GiNA(folder,resizeIndex,method,minApx,maxApx,realLength,keepOrder,showFig,writeT)
%
%
%       Cranberry Genetics and Genomics Lab, Horticulture Department
%                       University of Wisconsin-Madison
%
%                    if you use GINA, please cite:
%
% Diaz-Garcia L, Covarrubias-Pazaran G, Schlautman B, Zalapa J (2016) 
% GiNA, an Efficient and High-Throughput Software for Horticultural Phenotyping. 
% PLoS ONE 11(8): e0160439. https://doi.org/10.1371/journal.pone.0160439
%
disp(sprintf('\n\n\nThanks for using GiNA!\n\n\n'))

if nargin ~= 9
    error('GINA error. Number of argument inputs must be 9');
end

if isdir(folder)==0
    disp('A single image will be analyzed')
    data=miniGiNA(folder,'',resizeIndex,method,minApx,maxApx,realLength,keepOrder,showFig)
    data.file=folder;
    % writing excel
    if writeT
        WriteTGiNA(data, realLength)
    end
elseif isdir(folder)==1
    list=dir(fullfile(folder,'*.JPG'));
    tmpFolder=pwd;
    cd(folder)
    disp(strcat({'A folder containing '},num2str(length(list)),{' pictures will be analyzed'}))
    for superKK=1:length(list)
        
        folderX=list(superKK).name;
        folderY=list(superKK).folder
        
        data(superKK)=miniGiNA(folderX,folderY,resizeIndex,method,minApx,maxApx,realLength,keepOrder,showFig)
        
        disp(strcat(folderX,{' done'}))
    end
    % writing excel
    if writeT
        WriteTGiNA(data, realLength)
    end
end
end

function data=miniGiNA(folder,dir,resizeIndex,method,minApx,maxApx,realLength,keepOrder,showFig)
disp(strcat({'starting '},folder,'...'))
mark1=imresize(imread(folder),resizeIndex);
H = fspecial('disk', 4);
mark2=imfilter(mark1,H);

if isstruct(method)==1
    tmpTitle={'neural network implementation'};
    
    clear mark3
    m1=size(mark2);
    
    for i=1:3
        mark3(:,i)=reshape(mark2(:,:,i),m1(1)*m1(2),1);
    end
    
    a=method.network(double(mark3'));
    
    clear a3;clear mark5;clear mark4;
    a4=find(a(1,:)>0.4);   % set this value for improve background detection
    mark4=mark3;
    mark4(a4,:)=0;
    m1=size(mark2);
    
    for i=1:3
        mark5(:,:,i)=reshape(mark4(:,i),m1(1),m1(2));
    end
    
    mark6=mark5;
        
    [f,g]=find(mark5(:,:,2)>1);
    for i=1:3
        for j=1:length(f)
            mark6(f(j),g(j),i)=255;
        end
    end
    
    mark6=double(rgb2gray(mark6));
else
    tmpTitle=(strcat('threshold implementation, threshold set at',{' '},num2str(method(1))));
    mark2_beta=rgb2hsv(mark2);    %%%% check this later. exploring other color spaces
    mark6=double(mark2(:,:,method(2))<method(1));
    
end




% IMAGE ANALYSIS

clear B
clear binaryImage
sOnes=length(find(mark6==1));
sZeros=length(find(mark6==0));
if sOnes>sZeros
    B=1+(mark6*-1);
    warNeg='col conv implemented)';
else
    B=mark6;
    warNeg='col conv was not requiered)';
end

se = strel('disk',3);
binaryImage = imclose(B,se);
binaryImage = imclearborder(bwareaopen(binaryImage,minApx));
if showFig
    if isstruct(method)==0
        figure(1)
        subplot(3,1,1)
        imagesc(mark1)
        title('original resized picture')
        subplot(3,1,3)
        imagesc(binaryImage)
        title(strcat(tmpTitle,{', ('},warNeg))
        hold off
        subplot(3,1,2)
        imagesc(mark2(:,:,method(2)))
        title(strcat({'layer: '},num2str(method(2))))
        colorbar
        hold off
    end
    if isstruct(method)==1
        figure(1)
        subplot(3,1,1)
        imagesc(mark1)
        title('original resized picture')
        subplot(3,1,3)
        imagesc(binaryImage)
        title(strcat(tmpTitle,{', ('},warNeg))
        hold off
        subplot(3,1,2)
        imagesc(mark2(:,:,3))
        title(strcat({'layer: '},num2str(3)))
        colorbar
        hold off
    end
    [filepath,name,ext] = fileparts(folder);
    if ~strcmp(dir, '')
        filepath = dir
    end
    newfile = strcat(filepath,'/',name,'.binary.png');
    imwrite(binaryImage, newfile, 'BitDepth', 1);
    newfile = strcat(filepath,'/',name,'.segments',ext);
    mfile = strcat('export_fig',{' '},newfile,' -m2');
    eval(mfile{1});
end

labeledImage = bwlabel(binaryImage, 8);
blobMeasurements = regionprops(labeledImage, binaryImage, 'all');
numberOfBlobs = size(blobMeasurements, 1);
boundaries = bwboundaries(binaryImage);
numberOfBoundaries = size(boundaries);



for zz=1:length(blobMeasurements)
     tmpAA(zz)=blobMeasurements(zz).Area;
end

grandes=find(tmpAA>maxApx);

blobMeasurements(grandes)=[];
numberOfBlobs = size(blobMeasurements, 1);
boundaries(grandes)=[];
numberOfBoundaries = length(boundaries);



if isstruct(method)==1
    myRefs=[1:method.refs(1) numberOfBlobs-method.refs(2)+1:numberOfBlobs];
    numberOfRefs=sum(method.refs);
end
if isstruct(method)==0
    myRefs=[1:method(3) numberOfBlobs-method(4)+1:numberOfBlobs];
    numberOfRefs=sum(method(3:4));
end


if keepOrder
    for k = 1 : numberOfBlobs
        mygps(k,1:2)=blobMeasurements(k).BoundingBox(1:2);
    end
    [a,b]=sort(mygps(:,1),'ascend');
    
    
    idx=((numberOfRefs/2)+1):5:(numberOfBlobs-(numberOfRefs/2)-1);
    mygps2=[];
    for k=1:length(idx)
        tmpI=b(idx(k):idx(k)+4);
        [a2,b2]=sort(mygps(tmpI,2),'ascend');
        for l=1:length(a2)
            mygps2=[mygps2 tmpI(b2(l))];
        end
    end
    
    blobMeasurements2=blobMeasurements;
    boundaries2=boundaries;
    cont=1+numberOfRefs/2;
    for k=1:(numberOfBlobs-numberOfRefs)
        blobMeasurements2(cont)=blobMeasurements(mygps2(k));
        boundaries2{cont}=boundaries{mygps2(k)};
        cont=cont+1;
    end
else
    blobMeasurements2=blobMeasurements;
    boundaries2=boundaries;
    
end


blobMeasurements=blobMeasurements2;
boundaries=boundaries2;



for k = 1 : numberOfBlobs
    
    data.LvsW(k) = blobMeasurements(k).MajorAxisLength/blobMeasurements(k).MinorAxisLength;
    data.blobLength(k) = blobMeasurements(k).MajorAxisLength;
    data.blobWidth(k) = blobMeasurements(k).MinorAxisLength;
    data.projectedArea(k) = blobMeasurements(k).Area;
    data.projectedPerimeter(k) = blobMeasurements(k).Perimeter;
    B=0.5*blobMeasurements(k).MajorAxisLength;
    A=blobMeasurements(k).MinorAxisLength;
    C=B;
    data.skinSurface(k)=(2*pi*A^2)+(pi*A)*((((B^2)/sqrt((B^2)-(A^2)))*acosd(A/B))+(((C^2)/sqrt((C^2)-(A^2)))*acosd(A/C)));
    data.blobVolume(k)=((2*pi)/3)*(A^2)*(B+C);
    data.blobEccentricity(k)=blobMeasurements(k).Eccentricity;
    data.blobSolidity(k)=blobMeasurements(k).Solidity;
    data.locationX_All(k)=blobMeasurements(k).Centroid(1);
    data.locationY_All(k)=blobMeasurements(k).Centroid(2);
    
    
    clear tmp;
    
    for a=1:3
        for b=1:length(blobMeasurements(k).PixelList)
            tmp(1,b,a)=mark1(blobMeasurements(k).PixelList(b,2),blobMeasurements(k).PixelList(b,1),a);
        end
        
    end
    
    tmp900=double(rgb2gray(tmp));
    tmp901=median(tmp);
    tmp902=var(double(tmp));
    
    
    data.RGBcolor(k,:)=[tmp901(1,1,1) tmp901(1,1,2) tmp901(1,1,3) tmp902(1,1,1) tmp902(1,1,2) tmp902(1,1,3)];
    data.bwColor(k)=sum(1-tmp900/255)/data.projectedArea(k);
    data.vbwColor(k)=var(1-tmp900/255)*100;
    
end


fR=myRefs;
fF=setdiff(1:numberOfBlobs,fR);
data.numbering=[1:length(fF)];
data.LvsW_r=data.LvsW(fF)/median(data.LvsW(fR));
data.blobLength_r=data.blobLength(fF)/median(data.blobLength(fR));
data.blobWidth_r=data.blobWidth(fF)/median(data.blobWidth(fR));
data.projectedArea_r=data.projectedArea(fF)/median(data.projectedArea(fR));
data.projectedPerimeter_r=data.projectedPerimeter(fF)/median(data.projectedPerimeter(fR));
data.skinSurface_r=data.skinSurface(fF)/median(data.skinSurface(fR));
data.blobVolume_r=data.blobVolume(fF)/median(data.blobVolume(fR));
medEccentricity = median(data.blobEccentricity(fR));
data.blobEccentricity_r=data.blobEccentricity(fF)/median(data.blobEccentricity(fR));
data.blobSolidity_r=data.blobSolidity(fF)/median(data.blobSolidity(fR));
data.bwColor_r=data.bwColor(fF)/median(data.bwColor(fR));
data.vbwColor_r=data.vbwColor(fF)/median(data.vbwColor(fR));

%Reset these to only include the fruits (exclude references -- otherwise, the exported table output has the wrong rgb values)
data.RGBcolor = data.RGBcolor(fF,:)


% transformation px 2 in/cm
if islogical(realLength)==0
    digitalLength=median(data.blobLength(fR));
    ratioConv=realLength/digitalLength;
    fprintf('Size Reference Median Digital Length: %d pixels, Convert Ratio: %f cm/pixels\n', digitalLength, ratioConv);
    
    data.blobLength_r2=data.blobLength(fF)*ratioConv;
    data.blobWidth_r2=data.blobWidth(fF)*ratioConv;
    data.projectedArea_r2=data.projectedArea(fF)*(ratioConv^2);
    data.projectedPerimeter_r2=data.projectedPerimeter(fF)*ratioConv;
    data.skinSurface_r2=data.skinSurface(fF)*(ratioConv^2);
    data.blobVolume_r2=data.blobVolume(fF)*(ratioConv^3);
    data.accuracy=data.blobLength(fR)*ratioConv;
end
if islogical(realLength)==1
    data.accuracy=data.blobLength(fR)*ratioConv;
end


data.file=folder;
data.locationX=data.locationX_All(fF);
data.locationY=data.locationY_All(fF);

if showFig
    figure(2)
    imagesc(mark1);
    title(strcat(num2str(length(fF)),{' fruits included by using '},tmpTitle));
    hold on;
    
    
    for k = 1 : numberOfBlobs
        text(blobMeasurements(k).BoundingBox(1) , blobMeasurements(k).BoundingBox(2), num2str(k), 'FontSize', 14, 'FontWeight', 'Bold','color','r');
    end
    
    for k = 1 : numberOfBlobs
        thisBoundary = boundaries{k};
        plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 1);
    end
    
    for k=1:length(blobMeasurements)
        
        phi = linspace(0,2*pi,50);
        cosphi = cos(phi);
        sinphi = sin(phi);
        
        
        xbar = blobMeasurements(k).Centroid(1);
        ybar = blobMeasurements(k).Centroid(2);
        
        a = blobMeasurements(k).MajorAxisLength/2;
        b = blobMeasurements(k).MinorAxisLength/2;
        
        theta = pi*blobMeasurements(k).Orientation/180;
        R = [ cos(theta)   sin(theta)
            -sin(theta)   cos(theta)];
        
        xy = [a*cosphi; b*sinphi];
        xy = R*xy;
        
        x = xy(1,:) + xbar;
        y = xy(2,:) + ybar;
        
        
        
        line([x(1),x(25)],[y(1),y(25)],'LineWidth',2,'Color','m')
        line([x(13),x(38)],[y(13),y(38)],'LineWidth',2,'Color','c')
        %text(blobMeasurements(k).BoundingBox(1)+blobMeasurements(k).BoundingBox(3)+10,blobMeasurements(k).Centroid(2)-20,strcat('Area=',num2str(blobMeasurements(k).Area)),'FontSize', 10, 'color','b')
        %text(blobMeasurements(k).BoundingBox(1)+blobMeasurements(k).BoundingBox(3)+10,blobMeasurements(k).Centroid(2),strcat('PerimEter=',num2str(blobMeasurements(k).Perimeter)),'FontSize', 10, 'color','b')
        %text(blobMeasurements(k).BoundingBox(1)+blobMeasurements(k).BoundingBox(3)+10,blobMeasurements(k).Centroid(2)+20,strcat('BW color index=',num2str(data.bwColor(k))),'FontSize', 10, 'color','b')
        plot(blobMeasurements(k).Centroid(1),blobMeasurements(k).Centroid(2),'o','MarkerFaceColor','w','MarkerSize',3)
        hold on
    end
    [filepath,name,ext] = fileparts(folder);
    if ~strcmp(dir,'')
        filepath = dir
    end
    newfile = strcat(filepath,'/',name,'.blobs',ext);
    mfile = strcat({'export_fig '},newfile,{' -m2'});
    eval(mfile{1});
end
end

function FillMatrixGiNA(data, fields, fieldsTable)
    cont=1;
    clear filesNamesA; clear matrixA; clear locoK; clear tableStr;
    fieldsTableExpanded = {'''picture'''}
    for i = 1:length(fieldsTable)
        if iscell(fieldsTable{i})
            for j = 1:length(fieldsTable{i})
                fieldsTableExpanded{end+1} = strcat('''',fieldsTable{i}{j},'''')
            end
        else
            fieldsTableExpanded{end+1} = strcat('''',fieldsTable{i},'''')
        end
    end
    fieldsTableExpanded{end+1} = strcat('''','locationX','''')
    fieldsTableExpanded{end+1} = strcat('''','locationY','''')
    fieldsTableExpandedNew     = strcat('{',strjoin(fieldsTableExpanded, ','),'}')
    for ii=1:length(data)
        for i=1:length(data(ii).LvsW_r)
            cont_col=1;
            for j=1:length(fields)
                if strcmp(fields{j},'RGBcolor')
                    for k=1:length(fieldsTable{j})
                        matrixA(cont,cont_col)=data(ii).(fields{j})(i,k)
                        cont_col=cont_col+1
                    end
                else
                    matrixA(cont,cont_col)=data(ii).(fields{j})(i)
                    cont_col=cont_col+1
                end
                filesNamesA{cont,1}=data(ii).file
            end
            locoK(cont,1)=data(ii).locationX(i)
            locoK(cont,2)=data(ii).locationY(i)
            cont=cont+1
        end
    end
    tableStr = 'cFilt=table(filesNamesA,'
    [no_rows,no_cols] = size(matrixA)
    for i = 1:no_cols
        tableStr = strcat(tableStr,'matrixA(:,',num2str(i),'),')
    end
    tableStr =  strcat(tableStr, 'locoK(:,1),locoK(:,2),','''','VariableNames','''',',',fieldsTableExpandedNew, ')');
    eval(tableStr)
    writetable(cFilt,'GiNA_output.csv')
end

function WriteTGiNA(data, realLength)
    %What fields to export to the table, if writeT true.
    exportFields                  = {'numbering', 'LvsW_r','blobLength_r','blobWidth_r','projectedArea_r','projectedPerimeter_r','skinSurface_r','blobVolume_r','blobEccentricity_r','blobSolidity_r','bwColor_r','vbwColor_r','blobLength_r2','blobWidth_r2','projectedArea_r2','projectedPerimeter_r2','skinSurface_r2','blobVolume_r2','RGBcolor'}
    exportFields_realLength       = {'numbering', 'LvsW_r','blobLength_r','blobWidth_r','projectedArea_r','projectedPerimeter_r','skinSurface_r','blobVolume_r','blobEccentricity_r','blobSolidity_r','bwColor_r','vbwColor_r','RGBcolor'}
    exportFieldsTable             = {'numbering', 'LvsW_r','blobLength_r','blobWidth_r','projectedArea_r','projectedPerimeter_r','skinSurface_r','blobVolume_r','blobEccentricity_r','blobSolidity_r','bwColor_r','vbwColor_r','blobLength_r2','blobWidth_r2','projectedArea_r2','projectedPerimeter_r2','skinSurface_r2','blobVolume_r2',{'R_med','G_med','B_med','R_var','G_var','B_var'}}
    exportFieldsTable_realLength  = {'numbering', 'LvsW_r','blobLength_r','blobWidth_r','projectedArea_r','projectedPerimeter_r','skinSurface_r','blobVolume_r','blobEccentricity_r','blobSolidity_r','bwColor_r','vbwColor_r',{'R_med','G_med','B_med','R_var','G_var','B_var'}}
    if islogical(realLength)==0
        FillMatrixGiNA(data, exportFields, exportFieldsTable)
    elseif islogical(realLength)==1
        FillMatrixGiNA(data, exportFields_realLength, exportFieldsTable_realLength)
    end
end

