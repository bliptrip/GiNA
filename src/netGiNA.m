
function mynet=netGiNA(setBackground,varargin)

setBkg=imresize(imread(setBackground),.2);
m=size(setBkg);
myLayers=1+length(varargin);
for i=1:3
setBkg_reshape(:,i)=reshape(setBkg(:,:,i),m(1)*m(2),1);
end
setBkg_forCat=zeros(m(1)*m(2),myLayers);
setBkg_forCat(:,1)=1;
subplot(myLayers,1,1)
imagesc(setBkg)
title('background color for NN training')

superCat=[];
for i=1:length(varargin)
    subplot(myLayers,1,i+1)
    setFruits(i).image=imresize(imread(varargin{i}),.2);
    m=size(setFruits(i).image);
    for j=1:3
    setFruits(i).reshape(:,j)=reshape(setFruits(i).image(:,:,j),m(1)*m(2),1);
    end
    setFruits(i).forCat=zeros(m(1)*m(2),myLayers);
    setFruits(i).forCat(:,i+1)=1;
    
    if i==1
    superCat=cat(1,setBkg_forCat,setFruits(i).forCat);
    superCat2=cat(1,setBkg_reshape,setFruits(i).reshape);
    else
    superCat=cat(1,superCat,setFruits(i).forCat);
    superCat2=cat(1,superCat2,setFruits(i).reshape);
    end
    
    imagesc(setFruits(i).image)
    title(strcat('fruit color number',{' '},num2str(i),{' for NN training'}))
end


x = double(superCat2);
t = superCat;
trainFcn = 'trainscg';  

hiddenLayerSize = 10;
net = patternnet(hiddenLayerSize);

net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Train the Network
[net,tr] = train(net,x',t');
mynet.network=net;
