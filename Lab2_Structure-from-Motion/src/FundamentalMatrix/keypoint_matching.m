function [matches, fa, fb,da,db,score] = keypoint_matching(Ia, Ib, threshold)

peak_thresh = threshold;

Ia_s = single(Ia);
Ib_s = single(Ib);

[fa, da] = vl_sift(Ia_s, 'PeakThresh', peak_thresh) ;
[fb, db] = vl_sift(Ib_s, 'PeakThresh', peak_thresh) ;
[matches, score] = vl_ubcmatch(da, db) ;

end
