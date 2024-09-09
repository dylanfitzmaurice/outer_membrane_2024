function imout=norm16bit(imin,percent)
%written by Rico Rojas
%normalizes the intensity of a 16 bit image, saturating a certain 
%percentage of the pixels at the max and min values.

[imcounts,imbin]=hist(double(nonzeros(imin)),1000);
%figure,plot(imbin,imcounts),pause

csic=cumsum(fliplr(imcounts));
pcsic=csic/sum(imcounts);
[~,mpos1]=min(abs(pcsic-percent/100));
maxintensity=imbin(end-mpos1);

csic=cumsum(imcounts);
pcsic=csic/sum(imcounts');
[~,mpos2]=min(abs(pcsic-percent/100));

minintensity=imbin(mpos2);

imout=imadjust(imin,[minintensity/65535 maxintensity/65535],[]);