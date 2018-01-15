function blocksize = sft_setblocksize(patestimate,umperpix)

% Sets blocksize to include ~8 patterns.  This was found to be an optimal
% number of patterns within each scan to give an adequate, but localized,
% representation of patterning in the region.

disp(' ')
disp('Setting blocksize to include ~8 pattern widths within each frame')

blocksize = 2*(floor((patestimate*10/umperpix)/2));

% Make blocksize non-prime.
while isprime(blocksize)==1
    blocksize=blocksize+1;
end
disp(['Setting BLOCKSIZE = ' num2str(blocksize)]);

end