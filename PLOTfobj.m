function Z0=PLOTfobj( Chrom )
    % parameter estimation for 
    %Nind=pop_size;
    %Nvar=Dim
    load('skin_impedanceData.mat')
    [Nind,Nvar] = size(Chrom);
    Rd  = repmat(Chrom(:,1),1,size(f,2));
    Cd  = repmat(Chrom(:,2),1,size(f,2));
    Cp  = repmat(Chrom(:,3),1,size(f,2));
    a   = repmat(Chrom(:,4),1,size(f,2));
    b   = repmat(Chrom(:,5),1,size(f,2));
    c   = repmat(Chrom(:,6),1,size(f,2));
    r   = repmat(Chrom(:,7),1,size(f,2));
    
    one = ones(1,size(f,2));
%     Rp = a + b.*(one-exp(c.*f));
%     Rp = 1./(a + b.*(one-exp(c.*f)));
    Rp = a + b.*exp(c.*(f).^-1);
    
    Z_real = r+Rd./(one+(2*pi*f.*Rd.*Cd).^2)+Rp./(one+(2*pi*f.*Rp.*Cp).^2);
    Z_image =   (2*pi*f.*(Rd).^2.*Cd)./(one+(2*pi*f.*Rd.*Cd).^2)+...
                (2*pi*f.*(Rp).^2.*Cp)./(one+(2*pi*f.*Rp.*Cp).^2);
    Z0  = Z_real.^2 + Z_image.^2;
    Z0 = sqrt(Z0);
end