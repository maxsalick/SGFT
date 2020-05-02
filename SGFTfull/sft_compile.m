function [data, m_im, m_orig, m_cov, m_sl, m_dir, m_quiver, m_str, fillbins, superiorang, OI, AI, CMI,percsarc] = ...
    sft_compile(m_full_cov, imagesize, blocksize, m_full_sl, m_full_dir, ...
    m_full_str, quiver, im, imagedataorig, scanjump, umperpix)

% Takes output matrices and compiles them into data arrays.  Produces
% figures to report the resutls.

disp(' ')
disp('Producing report of results...')


m_cov = m_full_cov(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2));
ind = find(m_cov);
nonind = find(1-m_cov);




m_sl = m_full_sl(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2));
m_dir  = m_full_dir(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2));
m_quiver = quiver(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2));
m_str = (m_full_str(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2)).^1);
loweststrength = min(m_str(ind));
m_str = m_str + loweststrength.*(1-m_cov);

m_im = (im(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2)));
m_orig = (255-imagedataorig(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2)));


m_full_dist = umperpix.*(bwdist(1-m_full_cov));
m_dist = m_full_dist(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2));

% Make sparse matrices

m_sp_dist = m_dist(1:scanjump:end,1:scanjump:end);
m_sp_str = m_str(1:scanjump:end,1:scanjump:end);
m_sp_dir = m_dir(1:scanjump:end,1:scanjump:end);
m_sp_sl = m_sl(1:scanjump:end,1:scanjump:end);
m_sp_cov = m_cov(1:scanjump:end,1:scanjump:end);
[coveredY,coveredX] = size(m_sp_dist);
[Ymesh,Xmesh] = meshgrid(1:coveredY,1:coveredX);
Ymesh = Ymesh'; 
Ymesh = coveredY-Ymesh;
Xmesh = Xmesh';


lanewidth = (umperpix.*scanjump) .* max(max((m_sp_cov.*Xmesh)') - min((m_sp_cov.*Xmesh)'));
% laneedges = diff(m_full_cov)
[leng,wid]=size(m_full_cov);
for widthrow = 1:leng
    checkrow=[0 m_full_cov(widthrow,:) 0];
    rowends=strfind(checkrow,[1 0])-1;
	rowstrts=strfind(checkrow,[0 1]);
    if length(rowends)==1 && length(rowstrts)==1
    pixwidth(widthrow)=max(rowends-rowstrts+1);
    else
        pixwidth(widthrow)=0;
    end
end
lanewidth = mean(pixwidth(pixwidth>0))*umperpix;

data = [Xmesh(m_sp_cov>0.5) ...  % 1
                Ymesh(m_sp_cov>0.5)...  % 2
                m_sp_dist(m_sp_cov>0.5) ...  % 3
                m_sp_str(m_sp_cov>0.5) ...  %4
                m_sp_dir(m_sp_cov>0.5) ... %5
                m_sp_sl(m_sp_cov>0.5) ...  %6
                m_sp_cov(m_sp_cov>0.5) .* lanewidth];  %7


%% Alignment Distribution


[n,bins] = histc(data(:,5).*180./pi,[0:1:180]);

fillbins = zeros(size(n));

for k =1:max(size(data))
    fillbins(bins(k)) = fillbins(bins(k)) + (data(k,4)); 
end
fillbins(end)=[];

[angpks,anglocs] = findpeaks(fillbins);
[toss,highestangpk] = max(angpks);
allangX = sum(fillbins'.*cos(pi.*[2:2:360]./180));
allangY = sum(fillbins'.*sin(pi.*[2:2:360]./180));
if allangX>0
    superiorang = atan(allangY/allangX)/2;
else
    superiorang = (atan(allangY/allangX) + pi)/2;
end

if superiorang<0
    superiorang = superiorang+pi;
end

superiorang = superiorang*180/pi;

MRL = ((sum(cos(data(:,5))).^2+sum(sin(data(:,5))).^2).^.5)/length(data(:,5));


disp(' ')
disp('Data compiled and reported in "data" array.')

OI = mean(m_str(m_cov));
AI = MRL;
CMI = OI*AI;

data(:,8) = m_sp_cov(m_sp_cov>0.5) .* OI;
data(:,9) = m_sp_cov(m_sp_cov>0.5) .* AI;
data(:,10) = m_sp_cov(m_sp_cov>0.5) .* CMI;
percsarc = (nnz(m_str(m_cov)>.5)) / (nnz(m_cov)) * 100;
end