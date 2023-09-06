%% -------------ACTIVITY CALCULATION----------------

mask=load_untouch_nii('mask.nii'); 
pet=load_untouch_nii('pet.nii');

%load_untouch: load the nifti file, but do not apply any changes that
%are indicated in the header

% NOTE: You need to have a mask and a pet image both with the same size
if size(mask.img)==size(pet.img)
    fprintf('The pet and the mask have the same size')
end

%% COUNT VOXELS IN A MASK

%mask.img(:) returns the elements of the matrix in column order
mask_vox = sum(mask.img(:));

%Since the value of the voxel inside the mask is 1 and outside is 0
%(binary), the sum of all voxels gives the number of voxels

%% CREATE A VECTOR OF PET VOXELS BASED ON A MASK

%Returns a mask_vox-by-1 matrix of zeros; column matrix
vector_pet_vox = zeros(mask_vox,1);

[dim1,dim2,dim3] = size(pet.img); %Or size(mask.img)
i0=1;
for i = 1:dim1
    for j = 1:dim2
        for k = 1:dim3
            if mask.img(i,j,k) == 1
                vector_pet_vox(i0) = pet.img(i,j,k);
                i0=i0+1;
            end
        end
    end
end

%Transversing the mask, if the voxel is equal to 1, the value corresponding
%to that position is stored in a column vector with a dimension equal to
%the total number of voxels in the mask.
%The vector only contains the values corresponding to the volume of the pet
%defined by the organ segment.

%% VECTOR STATISTICAL INFORMATION in Bq/ml

%Mean
mean_pet_mask = mean(vector_pet_vox);

%Standard desviation
std_pet_mask=std(vector_pet_vox);

%Median
median_pet_mask=median(vector_pet_vox);

%Minimum
min_pet_mask=min(vector_pet_vox);

%Maximum
max_pet_mask=max(vector_pet_vox);

%Variance
var_pet_mask=var(vector_pet_vox);

%% SEGMENT VOLUME

%Since the segment consists of cubic voxels with a volume equal to 64 mm^3,
%the total volume is equal to: number of voxels * volume of each voxel
seg_volume_mm3=mask_vox*(4^3);
seg_volume_ml=seg_volume_mm3*0.001;

%% SEGMENT MEAN (AND MEDIAN) ACTIVITY IN Bq/ TOTAL ACTIVITY IN THE MASK (ORGAN) AT "THAT MOMENT" IN Bq

mean_activity_mask=mean_pet_mask*seg_volume_ml;
median_activity_mask=median_pet_mask*seg_volume_ml;

%% SEGMENT MEAN (AND MEDIAN) ACTIVITY IN MBq/ TOTAL ACTIVITY IN THE MASK (ORGAN) AT "THAT MOMENT" IN MBq

%Since each voxel has a Bq/ml value and the initial activity is in MBq, it
%is advisable to pass the units of activity to MBq

total_activity_mask=mean_activity_mask*0.000001;
total__median_activity_mask=median_activity_mask*0.000001;


%% PERCENTAGE OF ACTIVITY IN THE SEGMENT

A0=214.6; %administered activity em MBq

ratio_mean_act=(total_activity_mask/A0)*100;
ratio_median_act=(total__median_activity_mask/A0)*100;