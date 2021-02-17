function dets = get2DDeterminants(mats)
    if nargin==0
        N = 100;
        mats = randn(100,2,2);
    end
    dets = mats(:,1,1).*mats(:,2,2) - mats(:,1,2).*mats(:,2,1);
end