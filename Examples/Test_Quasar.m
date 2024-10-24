%%Test Quasar deconvolution
clear
load('quasar.mat');
f=(f-min(f(:)))/(max(f(:))-min(f(:)));

%build a gaussian kernel and blur the test image f
gauss=fspecial('gaussian',15,3);
fmir=repmat(f,3);
blurryf=conv2(fmir,gauss,'same');
blurryf=blurryf(size(f,1)+1:2*size(f,1),size(f,2)+1:2*size(f,2));

%perform the deconvolution
u=L1_L2_2D(blurryf,gauss,1e7,10,10,0);

%plot results
figure(1);
subplot(1,3,1);imshow(f,[]);
subplot(1,3,2);imshow(blurryf,[]);
subplot(1,3,3);imshow(u,[]);
