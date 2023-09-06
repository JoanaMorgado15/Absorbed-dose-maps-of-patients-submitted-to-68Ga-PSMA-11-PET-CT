
% CALCULATE DOSE

dose = zeros(dim1,dim2,dim3);   % dim1 x dim2 x dim3 are the dimensions of the PET map 
for i=dim:(dim1-(dim-1))        % dim is the dimension of the "s_kernel"
    for j=dim:(dim2-(dim-1))
        for k=dim:(dim3-(dim-1))
            for ii=-(dim-1):(dim-1)
                for jj=-(dim-1):(dim-1)
                    for kk=-(dim-1):(dim-1)
                        dose(i,j,k) = dose(i,j,k) + (s_kernel(ii+dim, jj+dim, kk+dim) .* pet.img(i+ii, j+jj, k+kk));
                    end
                end
            end
        end
    end
end
