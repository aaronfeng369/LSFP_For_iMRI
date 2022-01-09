function gamma=CoeffOper(op,alpha,beta)
% This subroutine implement alpha op beta = gamma, where op = '-' or 'l0'.
Level=length(alpha);
[nD,nD1]=size(alpha{1});

for ki=1:Level
    for ji=1:nD
        for jj=1:nD
            if op=='-'
                gamma{ki}{ji,jj}=alpha{ki}{ji,jj}-beta{ki}{ji,jj};
            elseif strcmp(op,'l0')  
                gamma{ki}{ji,jj}=(alpha{ki}{ji,jj}.^2-beta{ki}{ji,jj}>0).*alpha{ki}{ji,jj};
            end
        end
    end
end