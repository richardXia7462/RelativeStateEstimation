function errMark = GetErrorMark(markF, mark_est)
    nFeature = size(mark_est,2);
    errMark = zeros(nFeature,1);
    for i=1:nFeature
        errMark(i) = norm(markF(:,i)-mark_est(:,i));
    end
