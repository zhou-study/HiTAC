function [AC,PR,RE,F1] = AC_PE_RE_1(true_labels, cluster_labels)
%ACCURACY Compute clustering accuracy using the true and cluster labels and
%   return the value in 'score'.
%
%   Input  : true_labels    : N-by-1 vector containing true labels
%            cluster_labels : N-by-1 vector containing cluster labels
%
%   Output : score          : clustering accuracy

% Compute the confusion matrix 'cmat', where
%   col index is for true label (CAT),
%   row index is for cluster label (CLS).


n = length(true_labels);


%%%%%%%%%%%%%%%%%%%%%%%%%
max_t = max(max(true_labels));
max_c = max(max(cluster_labels));
if max_t < max_c
    temp = max_c;   
    for i = max_t:-1:1
        ifind = find(true_labels==i);
        true_labels(ifind) = temp;
        temp = temp - 1;
    end 
end
if max_c < max_t
    temp = max_t;   
    for i = max_c:-1:1
        ifind = find(cluster_labels==i);
        cluster_labels(ifind) = temp;
        temp = temp - 1;
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%


n = length(true_labels);
cat = spconvert([(1:n)' true_labels ones(n,1)]);
cls = spconvert([(1:n)' cluster_labels ones(n,1)]);
cls = cls';
% Confusion Matrix
cmat = full(cls * cat);
 
% Hungarian Algorithm
[match, cost] = hungarian(-cmat);

l = size(match,2); 
for i = 1 : l
    ind = find(match(:,i) == 1);    
    cmat([i,ind],:) = cmat([ind,i],:) ;
    match([i,ind],:) = match([ind,i],:) ;
end


sum_A = zeros(l,1);
sum_p = zeros(l,1);
sum_r = zeros(l,1);
for i = 1 : l
    sum_A(i) = cmat(i,i) ;
    sum_p(i) = cmat(i,i)/sum(cmat(i,:)) ;
    if sum(cmat(:,i)) == 0
        sum_r(i) = 0; % 避免除以 0
    else
        sum_r(i) = cmat(i,i) / sum(cmat(:,i));
    end
end

AC = sum(sum_A)/length(true_labels);
PR = mean(sum_p);
RE = mean(sum_r);
F1 = 2*PR*RE/(PR+RE);


