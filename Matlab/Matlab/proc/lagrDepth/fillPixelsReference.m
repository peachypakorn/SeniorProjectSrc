function final_labels = fillPixelsReference(Il, final_labels, gamma_c, gamma_d, r_median, maxDLagr)

[m,n] = size(final_labels);

occPix = zeros(m,n);
occPix(final_labels<-maxDLagr) = 1;

% Streak-based filling of invalidated pixels from left
fillVals = -ones(1,m)*maxDLagr;
final_labels_filled = final_labels;
for col=1:n
    curCol = final_labels(:,col);
    curCol(curCol<-maxDLagr) = fillVals(curCol<-maxDLagr);
    fillVals(curCol>=-maxDLagr) = curCol(curCol>=-maxDLagr);
    final_labels_filled(:,col) = curCol;
end

% Streak-based filling of invalidated pixels from right
fillVals = -ones(1,m)*maxDLagr;
final_labels_filled1 = final_labels;
for col=n:-1:1
    curCol = final_labels(:,col);
    curCol(curCol<-maxDLagr) = fillVals(curCol<-maxDLagr);
    fillVals(curCol>=-maxDLagr) = curCol(curCol>=-maxDLagr);
    final_labels_filled1(:,col) = curCol;
end

% Choose minimum disparity
final_labels = max(final_labels_filled,final_labels_filled1);

% Weighted median filtering on occluded regions (post processing)
final_labels_smoothed = weightedMedianMatlab(Il,final_labels,r_median,gamma_c,gamma_d);
final_labels(occPix==1) = final_labels_smoothed(occPix==1);

end

