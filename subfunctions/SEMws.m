function Y = SEMws(X,DIM_c,DIM_s,cFLAG)
%SEMws Within-Subject Error Variance
% For a matrix X, SEM = SEMws(X) returns the within-subject standard error
% of the mean, as defined by Cousineau, D. (2005). Confidence intervals in
% within-subject designs: A simpler solution to Loftus and Masson’s method.
% Tutorials in quantitative methods for psychology, 1(1), 42-45. The input
% matrix X must be a matrix with dimensions dim_c corresponding to the
% within subject conditions and dim_s corresponding to the individual
% subjects.
%
% The across subject variance is defined as the difference between
% the mean across all conditions X1_bar for each subject and the mean
% across all conditions and all subjects Xbar, calculated for each subject.
% This variance is then subtracted from the individual datapoints in X,
% resulting in a measure of error variance which is not affected by the
% across subject error variance. This error variance is then divided by the
% square root of the number of subjects in the usual way to calculate the
% Standard Error of the Mean (SEMws)
%
% Y = SEMws(X) calculates the across subject variance defined as
% X1_bar - Xbar where X1_bar is the mean across all rows and Xbar is the
% mean across all rows and columns. The within-subject SEM is then
% calculated as Y = std(X - X1_bar + Xbar)/sqrt(n) and returned as a column
% vector of size nx1 where n is the number of rows in X.
%
% Y = SEMws(X,DIM_c,DIM_s) uses the dimension DIM_c for the mean across
% conditions for each subject and dimenion DIM_s for the mean across
% subjects of the within-subject condition averages. The output is returned
% as a matrix with the same dimensions as X
%
% Y = SEMws(X,DIM_c,DIM_s,FLAG) will use the correction recommended by
% Morey (2008) which multiplies the standard error by M/(M-1), where M is
% the number of conditions. It has been shown that this correction brings
% the calculation of the within-subject SEM closer to the values produced
% by the original recommendation of Loftus and Masson (1994). The
% correction is performed if FLAG is set to 1, and not performed if FLAG is
% set to 0. The default if FLAG is omitted is 0.
%
% Example: If X = [1 2 3; 4 5 7; 6 9 10] where DIM_c corresponds to the
% rows of X and DIM_s corresponds to the columns of X
%
% then SEMws(X,1,2) = [0.2940 0.222 0.4006] and SEM(X,1,2,1) = [0.4410
% 0.3333 0.6009]
%
% Sayeed A.D. Kizuk March 2017

isDim1Set = nargin > 1 && ~ischar(DIM_c);
isDim2Set = nargin > 2 && ~ischar(DIM_s);
isFlagSet = nargin == 4 && ~ischar(FLAG);

DIM_X = numel(size(X));

correction = 1;

if ~(isDim1Set && isDim2Set)
    
    DIM_c = 1;
    DIM_s = 2;

if (isDim1Set || isDim2Set)
    
    error('Please specify both the condition dimension and the subject dimension');
    
end

end

if ~(DIM_c <= DIM_X && DIM_s <= DIM_X)

    error('The dimensions provided do not match the dimensions of X');
    
elseif isFlagSet == 1
    
    if ~or(cFLAG == 1,cFLAG == 0);
        
        error('The flag has to be either 1 or 0');
        
        if cFLAG == 1
            
            correction = size(X,DIM_c)/(size(X,DIM_c)-1);
            
        end
    end
else
    
    X1_bar = nanmean(X,DIM_c);
    Xbar = nanmean(nanmean(X,DIM_c),DIM_s);
    
    n = size(X,DIM_s);
    
    X_cen = bsxfun(@plus,bsxfun(@minus,X,X1_bar),Xbar);
    Y = nanstd(X_cen,[],DIM_s)/sqrt(n) * correction;
    
end
