function [fourstrengthval,sarclength,meandirect] = sft_scan3(section,orientsection,gradstr,blocksize,umperpix)

%  Scans through the entire image to detect patterning.  First conducts
%  gradient analysis to determine directionality of the image.  Next, 1-D
%  Fourier transforms are taken in the determined directions, and peak
%  detection of averaged Fourier transforms is conducted to assess the
%  amount of organization at the scanned areas.  Spacing between patterns
%  is determined spatially using intensity profiles across the image in the
%  known direction of patterning.

% gradstr(orientsection==0)=0;
% gradstr(orientsection==-pi)=0;
% gradstr(orientsection==-pi/2)=0;
% gradstr(orientsection==pi/2)=0;


C = sum(sum(gradstr.*cos(orientsection.*2))); % Center of mass X (circular form)
S = sum(sum(gradstr.*sin(orientsection.*2)));
meandirect = atan2(S,C);
if meandirect<0
    meandirect = meandirect+2*pi;
end
meandirect = meandirect/2;

% subplot(2,2,1)
% imagesc(orientsection)
% subplot(2,2,2)
% imagesc(gradstr)
% subplot(2,2,3)
% hist(orientsection(gradstr~=0),100)
% subplot(2,2,4)
% gscatter(orientsection(gradstr~=0),gradstr(gradstr~=0))
% pause(.0005)




npixels = floor(blocksize/1.1);
clear fourlina
fourlina = {zeros(1,npixels) zeros(1,npixels) zeros(1,npixels) zeros(1,npixels) zeros(1,npixels)};
warning('off','signal:findpeaks:noPeaks');
k=0;
fourlinsum = 0;
intensetot = 0;


for Doff=-8:1:8
    k=k+1;
    xi = [blocksize/2-blocksize*cos(meandirect)/2.2-Doff*sin(meandirect); ...
          blocksize/2+blocksize*cos(meandirect)/2.2-Doff*sin(meandirect)];
    yi = [blocksize/2+blocksize*sin(meandirect)/2.2-Doff*cos(meandirect); ...
          blocksize/2-blocksize*sin(meandirect)/2.2-Doff*cos(meandirect)];
      

       fourlina{k} = improfile(section,double(xi),double(yi),npixels);
       FOUR{k} = abs(fft(fourlina{k}));

       fourlinsum = fourlinsum+FOUR{k};
       intensetot = intensetot+sum(fourlina{k});
end



linesmade = k;

fourlin_ave = smooth(fourlinsum ./ linesmade);

fourlin_ave = fourlin_ave(1:ceil(0.5*npixels));
for k=length(fourlin_ave):-1:1
    fourlin_ave(k) = fourlin_ave(k)-min(fourlin_ave(1:k));
end

[pks,locs] = findpeaks(fourlin_ave,'sortstr','descend');

    if isempty(pks) == 0

%         Fdomain = locs(1);
    
        fourstrengthval = pks(1);
        
    else
        fourstrengthval = 0;

        
    end

str_four = fourstrengthval/intensetot;
if isnan(str_four) == 1
    str_four = 0;
end

%% 1D Profile Strength/Spacing Determination - Spacial Domain Analysis
% subaxis(2,2,2,1)

combodists = [];
fourlina2 = zeros(size(fourlina{1}));
for k=1:17
    fourlina2 = fourlina2+fourlina{k};
end
    fourlina2 = smooth(fourlina2);
%     plot(fourlina2)

%     hold on
    linemean = mean(fourlina2);
    linezeros = fourlina2 - linemean;
cross = find(diff(sign(linezeros)));

cross_l = linezeros(cross);
cross_r = linezeros(cross+1);

cross_spats = cross + cross_l./(cross_l-cross_r);
% scatter(cross_spats,ones(size(cross_spats)).*linemean)
%     pause(.00000001)
if length(cross_spats)>2
cross_dists = umperpix.*[diff(cross_spats(1:2:end)); diff(cross_spats(2:2:end))];
combodists = [combodists(:) ; cross_dists];

end
% end



if length(combodists)>1
    
combodists(combodists>mean(combodists)+1*std(combodists))=[];
combodists(combodists<mean(combodists)-1*std(combodists))=[];

sl_spatial = mean(combodists);
sl_spatial_std = std(combodists);
else
    sl_spatial = 0;
    sl_spatial_std = 0.1;
end
% hold off
% subaxis(2,2,2,2)

fourstrengthval = str_four*90;
sarclength = sl_spatial;
% sl_spatial;
% 
% bar([fourstrengthval-sl_spatial_std 1-abs(2-sl_spatial)])
% axis([0.5 2.5 0 4])



end