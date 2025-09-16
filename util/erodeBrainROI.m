function roi_brain_erode = erodeBrainROI( roi_brain, numPix, dilPix )
arguments
    roi_brain
    numPix = 3
    dilPix = 3
end
% Function to erode the ROI mask in the brain. numPix is the number of
% pixels to aim to erode at the end of the function. The function will
% first erode the ROI more than requested and then dilate back by the
% dilPix amount. This helps to create smoother ROI masks.

erPix = numPix+dilPix;
[iiEr,jjEr,kkEr] = ndgrid(-erPix:erPix, -erPix:erPix, -erPix:erPix);
nhoodEr = ( (iiEr.^2)/(erPix^2) + (jjEr.^2)/(erPix^2) + (kkEr.^2)/(erPix^2) <= 1);
SEer = strel(nhoodEr);
roi_brain_erode_pre_dil = imerode( logical(roi_brain), SEer );

if dilPix ~= 0
    [iiDil,jjDil,kkDil] = ndgrid(-dilPix:dilPix, -dilPix:dilPix, -dilPix:dilPix);
    nhoodDil = ( (iiDil.^2)/(dilPix^2) + (jjDil.^2)/(dilPix^2) + (kkDil.^2)/(dilPix^2) <= 1);
    SEdil = strel(nhoodDil);
    roi_brain_erode = imdilate( roi_brain_erode_pre_dil, SEdil );
else
    roi_brain_erode = roi_brain_erode_pre_dil;
end

end