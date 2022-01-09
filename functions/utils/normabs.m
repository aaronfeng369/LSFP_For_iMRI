function imgNorm = normabs (imgIn)
imgIn = abs(imgIn);
imgNorm = (imgIn - min(imgIn(:))) / (max(imgIn(:))-min(imgIn(:)));
end