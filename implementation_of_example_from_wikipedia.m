%generate data
S=[0.038 0.194 0.425 0.626 1.253 2.500 3.740];
Rate=[0.050 0.127 0.094 0.2122 0.2729 0.2665 0.3317];

beta=[0.9 0.2];i=1;
while i<6
    
    %calculate numerical outcome with beta
    for m=1:7
        numerical_values(m)=(beta(1)*S(m))/(beta(2)+S(m));
    end
    
    %calculate residue
    for j=1:7
        r(j)=Rate(j)-numerical_values(j);
    end
    rss(i)=sum(r.^2);
    r=transpose(r);
    while (size(r,2)==7)
        r=transpose(r);
    end
    
    %calculate jacobian
    diff1=(-S)./(beta(2)+S);
    diff2=(beta(1)*S)./((beta(2)+S).^2);
    change_of_r_1=transpose(diff1);
    change_of_r_2=transpose(diff2);
    J=[change_of_r_1 change_of_r_2];
    
    %put iteration here
    beta=transpose(beta);
    beta=beta-pinv(J)*r;
    beta=transpose(beta);
    i=i+1;
end