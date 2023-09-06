    clc; clear;
%clc: clear command window
%clear: remove items from workspace

%% LOAD FILES
pet_1=load_untouch_nii('CCC21339_PSMA_PET_1_60min.nii');
pet_2=load_untouch_nii('CCC21339_PSMA_PET_2_90min-transformed.nii');
liver=load_untouch_nii('Liver_Resample_CCC21339.nii');
kidneys=load_untouch_nii('Kidneys_Resample_CCC21339.nii');
spleen=load_untouch_nii('Spleen_Resample_CCC21339.nii');
%load_untouch: load the nifti file, but do not apply any changes that
%are indicated in the header.
%Uploading the pet images corresponding to the two acquisiton moments. The pet images have
%a Bq/ml concentration count.

%% VOXEL DIMENSIONS
pet_vox_dim = pet_1.hdr.dime.pixdim(1,2:4);
%Access to the header of the nifti pet file (each file has a header with
%metadata). In the header, access the "chapter" of the image dimensions and
%the float vector that contains the dimensions of the pixels. As it is a
%line vector, the first (and only) line is accessed and positions 2 to 4
%are stored in a vector, containing the voxel width, height and depth or
%slice thickness, respectively.

pet_vox_vol = 10^-3 * pet_vox_dim(1,1)*pet_vox_dim(1,2)*pet_vox_dim(1,3); 
%Then the volume of the voxel is calculated as the product of the three
%dimensions in cm (it is necessary to multiply to 10^-3), for the result to
%come in ml or cc (cubic centimeter) (the pet images have a Bq/ml concentration count).

%% CALIBRATE IMAGE
% remember that pet image voxels are expressed in Bq/ml

pet_1_MBq = pet_1.img * 2^(-30/68) * 10^-6 * pet_vox_vol;  
pet_2_MBq = pet_2.img * 2^(-75/68)* 10^-6 * pet_vox_vol;
% now pet voxels are in MBq
%Since pet images voxels are expressed in Bq/ml, it is necessary to
%multiply all voxels of the pet image by the volume to obtain the activity
%in each voxel in Bq. Since the initial activity is in Mbq, it is better to
%multiply by 10^-6 to obtain the units of activity in MBq.

%% S_kernel
% Read in s-values in data

s_values_tissue = importdata('Ga68_4p0_SoftTissue.txt','\t',2); 

%importdata loads data into array s_values, tab (\t) is the column
%separator in the ASCII file and reading numeric data starts from line 3
%(headerlineIn+1; 2+1=3).

s_values_tissue =s_values_tissue.data(:,4);  
%s_values_red=s_values_red.data(:,4);
% S-values are in the 4th collumn of the imported file, so they are saved
% in a column vector.

dim = nthroot(length(s_values_tissue),3); %Equal
% Returns de cubic root of the previews array. The total number of s_values
% can be seen as the octant volume. Thus, to know the dimensions of the
% octant, that is, width, length and depth, it is necessary to calculate
% the cubic root of the total number.
s_values_tissue=deal(zeros(dim, dim, dim));
%Creation of a three-dimensional matrix of zeros with dimension 4x4x4 in 
%order to be filled.

for n = 1:length(s_values_tissue)
    k = fix((n-1)/(dim*dim)); %Posição anterior na área do octante; de 0 a 3
    r = (n-1)-k*dim*dim; %Resto
    j = fix(r/dim); %O que resta num eixo; de 0 a 3
    i = r-j*dim; %O que sobra; r-r'; de 0 a 3
    
    j = j+1;
    k = k+1;
    i = i+1;
    s_values_tissue_img(i,j,k) = s_values_tissue(n);
end
%The objective is to go through all the positions of the s_values vector
%and store the values in an octant in the corresponding (i,j,k) positions.
%The fix function is used to round the numbers towards zero since the
%coordinates are integers. It also starts at n-1 because the coordinates
%start at 0, but the positions in the vector start at 1.
%Starts by filling the k, then j and i taking into account the "rest", what
%"remains", since fix rounds to the smallest value.


% Constructing the total 3D kernel

s_values_all_tissue=deal(zeros(2*dim-1, 2*dim-1, 2*dim-1));
%In order to move from an octant to a 3D kernel, a matrix is created
%with twice the coordinates/positions but it is subtracted 1 so that the
%central voxel/coordinate is not repeated. 4x4-1=7; -3 -2 -1 0 1 2 3; an
%octant 4x4x4 is going to become a 3D kernel with dimensions 7x7x7.


for i = 1:length(s_values_all_tissue) %Equal
    ii = abs(i-dim)+1;
    for j = 1:length(s_values_all_tissue)
        jj = abs(j-dim)+1;
        for k = 1:length(s_values_all_tissue)
            kk = abs(k-dim)+1;
            s_values_all_tissue(i,j,k) = s_values_tissue_img(ii,jj,kk);
        end
    end
end
%Since the same radiopharmaceutical is used in the same tissue (soft
%tissue), the concept of symmetry is used in order to fill the other 7
%octats from the first one. i,j,k=1:7 (7x7x7 kernel) and ii,jj,kk=1:4
%(4x4x4), abs returns the absolute volume. Then, the s-value calculated
%from the first octant in a given position is mirrored for the positions in
%the other octants and the 3D kernel consists of the s-values that will be
%used in the dose calculation. 

%% ABSORBED DOSE

% PET MAP DIMENSIONS: dim1 x dim2 x dim3

[dim1,dim2,dim3] = size(pet_1_MBq); %Or size(pet_2_MBq), the images have the same size

[dose_1_tissue, dose_2_tissue] = deal(zeros(dim1,dim2,dim3)); 
%Returns an dim1xdim2xdim3 array of zeros in order to be later filled

for i=dim:(dim1-(dim-1)) % 'dim' is each dimension of "S-dose convolution kernel" 
    for j=dim:(dim2-(dim-1))
        for k=dim:(dim3-(dim-1))
            for ii=-(dim-1):(dim-1)
                for jj=-(dim-1):(dim-1)
                    for kk=-(dim-1):(dim-1)
                        dose_1_tissue(i,j,k) = dose_1_tissue(i,j,k) + ( s_values_all_tissue(ii+dim, jj+dim, kk+dim) .* pet_1_MBq(i+ii, j+jj, k+kk));
                        dose_2_tissue(i,j,k) = dose_2_tissue(i,j,k) + (s_values_all_tissue(ii+dim, jj+dim, kk+dim) .* pet_2_MBq(i+ii, j+jj, k+kk));

                    end
                end
            end
        end
    end
end

%(i,j,k) are the coordinates of the pet voxels in which dose will be
%calculated. 
%(dim-1) due to the effect of peripheries. For each voxel to be calculated
%there are at least 3 voxels in each direction i,j,k. Thus, it's not
%possible to calculate the dose in the 3 most peripheral voxels, where
%there is a frame at 0. Only when we reach the 4th voxel there is enough
%information and from the 166th voxel this information no longer exists.
%The dose formula is present here, where we have the
%activity in the voxel-source and the s-values (contribution of the
%voxel-source to the voxel-target). These are multiplied since they are in
%the same position in space and are both 4mm. There is a summation because
%the dose in the voxel-target needs the contribution of each voxel-source
%at a distance considered depending on the size of the kernel that is
%stored in a variable. As the scan takes place, contributions are saved and
%added to the existing variable. 

%% COMPUTE THE TIME INTEGRAL

halfLife = 67.71 * 60; % Half life of the radionuclide expressed in seconds
t0 = 30 * 60;
t1 = 75 * 60;
t2 = Inf;

Ga68_integral_1 = calculateIntegral(halfLife,t0,t1);
Ga68_integral_2 = calculateIntegral(halfLife,t1,t2);

%% ABSOLUTE ABSORBED DOSE(Gy)

dose_Gy_1_tissue=deal(pet_1); 
dose_Gy_2_tissue=deal(pet_2); 

dose_Gy_1_tissue.img = dose_1_tissue * 10^-3 * Ga68_integral_1;
dose_Gy_2_tissue.img = dose_2_tissue * 10^-3 * Ga68_integral_2;


dose_Gy_total_tissue=deal(pet_1);
dose_Gy_total_tissue.img=dose_Gy_1_tissue.img+dose_Gy_2_tissue.img;

%We want the dose structure to have the same parameters as the initial
%image, including the header. However, the content of the image is changed
%to obtain the dose map that results from multiplying the dose calculated
%by the integral (definition of the decay time limits to obtain acumulated
%activity). Since we want the map in Gy and have it in mGy, is necessary to
%multiply by 10^-3.

save_untouch_nii(dose_Gy_total_tissue, 'pet_tissue_dose_4p0.nii')
% Save the structure dose_Gy with the dose map (but do not apply any changes that
% are indicated in the header) in the current folder with the name
% 'pet_dose_4p0.nii'.

%% MASK DOSE (Gy)


[liver_dose_1, heart_dose_1, kidneys_dose_1, lungs_dose_1, brain_dose_1]=deal(dose_Gy_1);
[liver_dose_2, heart_dose_2, kidneys_dose_2, lungs_dose_2, brain_dose_2]=deal(dose_Gy_2);
%We want the mask dose structure to have the same parameters as the total
%dose map,including the header and the dose (in Gy) in the image field.
%deal assings all structures at once.

[liver_vox_dose_1,liver_vox_dose_2] = deal(zeros(sum(liver.img(:)),1));
[heart_vox_dose_1, heart_vox_dose_2] = deal(zeros(sum(heart.img(:)),1));
[kidneys_vox_dose_1,kidneys_vox_dose_2] = deal(zeros(sum(kidneys.img(:)),1));
[lungs_vox_dose_1,lungs_vox_dose_2] = deal(zeros(sum(lungs.img(:)),1));
[brain_vox_dose_1,brain_vox_dose_2]  = deal(zeros(sum(brain.img(:)),1));
%mask.img(:) returns the elements of the matrix in column order
%Since the value of the voxel inside the mask is 1 and outside is 0
%(binary), the sum of all voxels gives the number of voxels
%Returns a mask_vox-by-1 matrix of zeros; column matrix

[dim1,dim2,dim3] = size(dose_Gy_1.img);
[i0, j0, k0, e0, r0]=deal(1);
for i = 1:dim1
    for j = 1:dim2
        for k = 1:dim3
            if liver.img(i,j,k) == 1
                liver_vox_dose_1(i0)=dose_Gy_1.img(i,j,k);
                liver_vox_dose_2(i0)=dose_Gy_2.img(i,j,k);
                i0=i0+1;
            else
                [liver_dose_1.img(i,j,k),liver_dose_2.img(i,j,k)]=deal(0);
            end
                if heart.img(i,j,k)==1
                    heart_vox_dose_1(j0)=dose_Gy_1.img(i,j,k);
                    heart_vox_dose_2(j0)=dose_Gy_2.img(i,j,k);
                    j0=j0+1;
                else
                    [heart_dose_1.img(i,j,k),heart_dose_2.img(i,j,k)]=deal(0);

                end
                    if kidneys.img(i,j,k)==1
                        kidneys_vox_dose_1(k0)=dose_Gy_1.img(i,j,k);
                        kidneys_vox_dose_2(k0)=dose_Gy_2.img(i,j,k);
                        k0=k0+1;
                    else
                        [kidneys_dose_1.img(i,j,k),kidneys_dose_2.img(i,j,k)]=deal(0);
                    end
                        if lungs.img(i,j,k)==1
                            lungs_vox_dose_1(e0)=dose_Gy_1.img(i,j,k);
                            lungs_vox_dose_2(e0)=dose_Gy_2.img(i,j,k);
                            e0=e0+1;
                        else
                            [lungs_dose_1.img(i,j,k),lungs_dose_2.img(i,j,k)]=deal(0);
                        end
                            if brain.img(i,j,k)==1
                                brain_vox_dose_1(r0)=dose_Gy_1.img(i,j,k);
                                brain_vox_dose_2(r0)=dose_Gy_2.img(i,j,k);
                                r0=r0+1;
                            else
                                [brain_dose_1.img(i,j,k),brain_dose_2.img(i,j,k)]=deal(0);
                            end             
        end
    end
end

%Transversing the mask, if the voxel is equal to 1, the value corresponding
%to that position is stored in a column vector with a dimension equal to
%the total number of voxels in the mask.
%The vector only contains the values corresponding to the volume of the pet
%defined by the organ segment.
%Since mask_dose.img is composed of the total body dose values, a condition is required
%so that outside the mask (when it's equal to 0), the values of
%mask_dose.img are equal to 0 (we only want the dose inside the mask)

[liver_dose_total, heart_dose_total, kidneys_dose_total, lungs_dose_total, brain_dose_total]=deal(dose_Gy_1);
liver_dose_total.img=liver_dose_1.img+liver_dose_2.img;
heart_dose_total.img=heart_dose_1.img+heart_dose_2.img;
kidneys_dose_total.img=kidneys_dose_1.img+kidneys_dose_2.img;
lungs_dose_total.img=lungs_dose_1.img+lungs_dose_2.img;
brain_dose_total.img=brain_dose_1.img+brain_dose_2.img;

save_untouch_nii(liver_dose_1, 'liver_dose_1_4p0.nii')
save_untouch_nii(liver_dose_2, 'liver_dose_2_4p0.nii')
save_untouch_nii(liver_dose_total, 'liver_dose_total_4p0.nii')
save_untouch_nii(heart_dose_1, 'heart_1_dose_4p0.nii')
save_untouch_nii(heart_dose_2, 'heart_2_dose_4p0.nii')
save_untouch_nii(heart_dose_total, 'heart_dose_total_4p0.nii')
save_untouch_nii(kidneys_dose_1, 'kidneys_dose_1_4p0.nii')
save_untouch_nii(kidneys_dose_2, 'kidneys_dose_2_4p0.nii')
save_untouch_nii(kidneys_dose_total, 'kidneys_dose_total_4p0.nii')
save_untouch_nii(lungs_dose_1, 'lungs_dose_1_4p0.nii')
save_untouch_nii(lungs_dose_2, 'lungs_dose_2_4p0.nii')
save_untouch_nii(lungs_dose_total, 'lungs_dose_total_4p0.nii')
save_untouch_nii(brain_dose_1, 'brain_dose_1_4p0.nii')
save_untouch_nii(brain_dose_2, 'brain_dose_2_4p0.nii')
save_untouch_nii(brain_dose_total, 'brain_dose_total_4p0.nii')

%Save the structure mask_dose with the dose map (but do not apply any changes that
%are indicated in the header) in the current folder with the name
%'mask_dose_4p0.nii'. 

%% MASK DOSE STATISTICAL INFORMATION (Gy)

organ={'Liver'; 'Heart'; 'Kidneys'; 'Lungs'; 'Brain'};
mean_1=[mean(liver_vox_dose_1)*1000; mean(heart_vox_dose_1)*1000; mean(kidneys_vox_dose_1)*1000; mean(lungs_vox_dose_1)*1000; mean(brain_vox_dose_1)*1000];
median_1=[median(liver_vox_dose_1)*1000; median(heart_vox_dose_1)*1000; median(kidneys_vox_dose_1)*1000; median(lungs_vox_dose_1)*1000; median(brain_vox_dose_1)*1000];
min_1=[min(liver_vox_dose_1)*1000; min(heart_vox_dose_1)*1000; min(kidneys_vox_dose_1)*1000; min(lungs_vox_dose_1)*1000; min(brain_vox_dose_1)*1000];
max_1=[max(liver_vox_dose_1)*1000; max(heart_vox_dose_1)*1000; max(kidneys_vox_dose_1)*1000; max(lungs_vox_dose_1)*1000; max(brain_vox_dose_1)*1000];
std_1=[std(liver_vox_dose_1)*1000; std(heart_vox_dose_1)*1000; std(kidneys_vox_dose_1)*1000; std(lungs_vox_dose_1)*1000; std(brain_vox_dose_1)*1000];
var_1=[var(liver_vox_dose_1)*1000; var(heart_vox_dose_1)*1000; var(kidneys_vox_dose_1)*1000; var(lungs_vox_dose_1)*1000;var(brain_vox_dose_1)*1000];

table_1=table(organ, mean_1, median_1, min_1, max_1, std_1, var_1);

mean_2= [mean(liver_vox_dose_2)*1000; mean(heart_vox_dose_2)*1000; mean(kidneys_vox_dose_2)*1000; mean(lungs_vox_dose_2)*1000; mean(brain_vox_dose_2)*1000];
median_2=[median(liver_vox_dose_2)*1000; median(heart_vox_dose_2)*1000; median(kidneys_vox_dose_2)*1000; median(lungs_vox_dose_2)*1000; median(brain_vox_dose_2)*1000];
min_2=[min(liver_vox_dose_2)*1000; min(heart_vox_dose_2)*1000; min(kidneys_vox_dose_2)*1000; min(lungs_vox_dose_2)*1000; min(brain_vox_dose_2)*1000];
max_2= [max(liver_vox_dose_2)*1000; max(heart_vox_dose_2)*1000; max(kidneys_vox_dose_2)*1000; max(lungs_vox_dose_2)*1000; max(brain_vox_dose_2)*1000];
std_2=[std(liver_vox_dose_2)*1000; std(heart_vox_dose_2)*1000; std(kidneys_vox_dose_2)*1000; std(lungs_vox_dose_2)*1000; std(brain_vox_dose_2)*1000];
var_2=[var(liver_vox_dose_2)*1000; var(heart_vox_dose_2)*1000; var(kidneys_vox_dose_2)*1000; var(lungs_vox_dose_2)*1000; var(brain_vox_dose_2)*1000];

table_2=table(organ, mean_2, median_2, min_2, max_2, std_2, var_2);

%Statistical analysis is easier with the mean, median, min, max, std and
%var commands. These have input arrays, that's why the mask dose values were
%preivously stored in a vector.
%Since the values are in the order of mGy, for a better visualization the
%conversion from Gy to mGy was made by multiplying by 1000.
%For a better visualization of the results, tables were created with the
%table command (table array with named variables that contain different
%types), since it was intended to have the numerical data with statistical
%values but also captions. 

mean_total=mean_1+mean_2;
median_total=median_1+median_2;
min_total=min_1+min_2;
max_total=max_1+max_2;
std_total=std_1+std_2;
var_total=var_1+var_2;
%What is intended is the statistical data (mainly mean and median) total,
%so they were determined separately in each distribution (pet1 and pet2)
%and then added for each segment (organ).

table_total=table(organ, mean_total, median_total, min_total, max_total, std_total, var_total);

writetable(table_1, 'dose_organ.xlsx', 'Sheet', 1);
writetable(table_2, 'dose_organ.xlsx', 'Sheet', 2);
writetable(table_total, 'dose_organ.xlsx', 'Sheet', 3);
%If we want to work with the data outside the Matlab, the tables are saved
%in an excel file with the name 'dose_organ.xlsx', each on a different
%sheet.

%% DOSE VOLUME HISTOGRAM

DVH=csvread('Segments_DVH.csv',1,0); %Reads a comma-separated value (CSV) formatted file (only numeric values) into array DVH. Reads data from
%the file starting at row 1 (the line 0 are titles).

%Since the volume data from Slicer are in %, to obtain the volume values in
%ml just multiply by the total volume (corresponding to 100%) and divide by
%100.
v100_kidneys_ml=316.032;v100_lungs_ml=2395.540;v100_liver_ml=1121.250;v100_heart_ml=887.208;v_100_brain_ml=1221.330;


kidneys_DVH_frac=[DVH(:,1)*1000 DVH(:,2)];
kidneys_DVH_plot_cum_frac=plot(kidneys_DVH_frac(:,1), kidneys_DVH_frac(:,2));%Creates a 2D line plot. DVH(:,1) are the elements of column 1 and DVH(:,2) of column 2
xlabel('Dose (mGy)');
ylabel('Volume Fraction (%)');
title('Cumulative DVH Kidneys Fraction');
saveas(kidneys_DVH_plot_cum_frac, 'Cumulative DVH Kidneys Fraction.png');
figure;
kidneys_DVH_ml=[DVH(:,1)*1000 (DVH(:,2)*v100_kidneys_ml)/100];
kidneys_DVH_plot_cum_ml=plot(kidneys_DVH_ml(:,1), kidneys_DVH_ml(:,2));
ylim([0 v100_kidneys_ml]);
xlabel('Dose (mGy)');
ylabel('Volume (ml)');
title('Cumulative DVH Kidneys Volume');
saveas(kidneys_DVH_plot_cum_ml, 'Cumulative DVH Kidneys Volume.png');
figure;
lungs_DVH_frac=[DVH(:,3)*1000 DVH(:,4)];
lungs_DVH_plot_cum_frac=plot(lungs_DVH_frac(:,1), lungs_DVH_frac(:,2));
xlabel('Dose (mGy)');
ylabel('Volume Fraction (%)');
title('Cumulative DVH Lungs Fraction');
saveas(lungs_DVH_plot_cum_frac, 'Cumulative DVH Lungs Fraction.png');
figure;
lungs_DVH_ml=[DVH(:,3)*1000 (DVH(:,4)*v100_lungs_ml)/100];
lungs_DVH_plot_cum_ml=plot(lungs_DVH_ml(:,1), lungs_DVH_ml(:,2));
ylim([0 v100_lungs_ml]);
xlabel('Dose (mGy)');
ylabel('Volume (ml)');
title('Cumulative DVH Lungs Volume');
saveas(lungs_DVH_plot_cum_ml, 'Cumulative DVH Lungs Volume.png');
figure;
liver_DVH_frac=[DVH(:,5)*1000 DVH(:,6)];
liver_DVH_plot_cum_frac=plot(liver_DVH_frac(:,1), liver_DVH_frac(:,2));
xlabel('Dose (mGy)');
ylabel('Volume Fraction (%)');
title('Cumulative DVH Liver Fraction');
saveas(liver_DVH_plot_cum_frac, 'Cumulative DVH Liver Fraction.png');
figure;
liver_DVH_ml=[DVH(:,5)*1000 (DVH(:,6)*v100_liver_ml)/100];
liver_DVH_plot_cum_ml=plot(liver_DVH_ml(:,1), liver_DVH_ml(:,2));
ylim([0 v100_liver_ml]);
xlabel('Dose (mGy)');
ylabel('Volume (ml)');
title('Cumulative DVH Liver Volume');
saveas(liver_DVH_plot_cum_ml, 'Cumulative DVH Liver Volume.png');
figure;
heart_DVH_frac=[DVH(:,7)*1000 DVH(:,8)];
heart_DVH_plot_cum_frac=plot(heart_DVH_frac(:,1), heart_DVH_frac(:,2));
xlabel('Dose (mGy)');
ylabel('Volume Fraction (%)');
title('Cumulative DVH Heart Fraction');
saveas(heart_DVH_plot_cum_frac, 'Cumulative DVH Heart Fraction.png');
figure;
heart_DVH_ml=[DVH(:,7)*1000 (DVH(:,8)*v100_heart_ml)/100];
heart_DVH_plot_cum_ml=plot(heart_DVH_ml(:,1), heart_DVH_ml(:,2));
ylim([0 v100_heart_ml]);
xlabel('Dose (mGy)');
ylabel('Volume (ml)');
title('Cumulative DVH Heart Volume');
saveas(heart_DVH_plot_cum_ml, 'Cumulative DVH Heart Volume.png');
figure;
brain_DVH_frac=[DVH(2:size(DVH,1),9)*1000 DVH(2:size(DVH,1),10)];
brain_DVH_plot_cum_frac=plot(brain_DVH_frac(:,1), brain_DVH_frac(:,2));
xlabel('Dose (mGy)');
ylabel('Volume Fraction (%)');
title('Cumulative DVH Brain Fraction');
saveas(brain_DVH_plot_cum_frac, 'Cumulative DVH Brain Fraction.png');
figure;
brain_DVH_ml=[DVH(2:size(DVH,1),9)*1000 (DVH(2:size(DVH,1),10)*v_100_brain_ml)/100];
brain_DVH_plot_cum_ml=plot(brain_DVH_ml(:,1), brain_DVH_ml(:,2));
ylim([0 v_100_brain_ml]);
xlabel('Dose (mGy)');
ylabel('Volume (ml)');
title('Cumulative DVH Brain Volume');
saveas(brain_DVH_plot_cum_ml, 'Cumulative DVH Brain Volume.png');

figure;
all_DVH_cum_frac=plot(DVH(:,1)*1000, DVH(:,2));
xlabel('Dose (mGy)');
ylabel('Volume Fraction (%)');
title('Cumulative DVH Segments Fraction');
hold on
plot(DVH(:,3)*1000, DVH(:,4));
hold on
plot(DVH(:,5)*1000, DVH(:,6));
hold on
plot(DVH(:,7)*1000, DVH(:,8));
hold on
plot(DVH(2:size(DVH,1),9)*1000, DVH(2:size(DVH,1),10));
hold off
legend('Kidneys', 'Lungs', 'Liver', 'Heart', 'Brain');
saveas(all_DVH_cum_frac, 'Cumulative DVH Segments Fraction.png');

mean_cDVH=[trapz(liver_DVH_frac(:,1), liver_DVH_frac(:,2))/100; trapz(heart_DVH_frac(:,1), heart_DVH_frac(:,2))/100; trapz(kidneys_DVH_frac(:,1), kidneys_DVH_frac(:,2))/100; trapz(lungs_DVH_frac(:,1), lungs_DVH_frac(:,2))/100; trapz(brain_DVH_frac(:,1), brain_DVH_frac(:,2))/100];
%The area under the cumulative histogram curve is equal to the mean dose of
%the segments. trapz(X,Y) integrates Y with respect to the coordinates or
%scalar spacing specified by X. Since the volume is in %, it is divided by
%100 to obtain the normalized value in 1.

% [kidneys_median_cDVH, kidneys_min_cDVH, lungs_median_cDVH,lungs_min_cDVH, liver_median_cDVH,liver_min_cDVH, heart_median_cDVH,heart_min_cDVH, brain_median_cDVH, brain_min_cDVH]=deal(zeros);%scalar zero
% 
% for i=1:length(kidneys_DVH_frac)
%     if kidneys_DVH_frac(i,2)>50 &&  kidneys_DVH_frac(i,2)<51
%         kidneys_median_cDVH=kidneys_DVH_frac(i,1);
%     end
%     if kidneys_DVH_frac(i,2)<100 && kidneys_DVH_frac(i,2)>98
%        kidneys_min_cDVH=kidneys_DVH_frac(i,1);
%     end
%     
%     if lungs_DVH_frac(i,2)>53 && lungs_DVH_frac(i,2)<54
%         lungs_median_cDVH=lungs_DVH_frac(i,1);
%     end
%     if lungs_DVH_frac(i,2)<100 && lungs_DVH_frac(i,2)>99
%         lungs_min_cDVH=lungs_DVH_frac(i,1);
%     end
%     if liver_DVH_frac(i,2)>49 && liver_DVH_frac(i,2)<50
%         liver_median_cDVH=liver_DVH_frac(i,1);
%     end
%     if liver_DVH_frac(i,2)<100 && liver_DVH_frac(i,2)>99.998
%         liver_min_cDVH=liver_DVH_frac(i,1);
%     end
%     if heart_DVH_frac(i,2)>50 && heart_DVH_frac(i,2)<51
%         heart_median_cDVH=heart_DVH_frac(i,1);
%     end
%     if heart_DVH_frac(i,2)<100 && heart_DVH_frac(i,2)>99.97
%         heart_min_cDVH=heart_DVH_frac(i,1);
%     end
% end
% for i=1:length(brain_DVH_frac)
%     if brain_DVH_frac(i,2)>51 && brain_DVH_frac(i,2)<52
%         brain_median_cDVH=brain_DVH_frac(i,1);
%     end
%     if brain_DVH_frac(i,2)<100 && brain_DVH_frac(i,2)>99.996
%         brain_min_cDVH(i)=brain_DVH_frac(i,1);
%         brain_min_cDVH_f=brain_min_cDVH(2);
%         
%     end
% end
% 
% median_cDVH=[liver_median_cDVH; heart_median_cDVH; kidneys_median_cDVH; lungs_median_cDVH; brain_median_cDVH];
% %For the median value, since the cDVH values do not provide a value of %
% %equal to 50%, approximations were made in this case and for each segment
% %individually. An approximation is also made in the min value.
% %When there is more than one option, a vector is created and at the end the
% %element is chosen
% 
% min_cDVH=[liver_min_cDVH; heart_min_cDVH; kidneys_min_cDVH; lungs_min_cDVH; brain_min_cDVH_f];
% max_cDVH=[max(liver_DVH_frac(:,1)); max(heart_DVH_frac(:,1)); max(kidneys_DVH_frac(:,1)); max(lungs_DVH_frac(:,1)); max(brain_DVH_frac(:,1))];
% 

table_cDVH=table(organ, mean_cDVH);
writetable(table_cDVH, 'dose_organ_DVH.xlsx', 'Sheet', 1);

figure;
all_DVH_cum_ml=plot(DVH(:,1)*1000, (DVH(:,2)*v100_kidneys_ml)/100);
ylim([0 v100_lungs_ml]); %Since the volume of the lungs is the largest volume
xlabel('Dose (mGy)');
title('Cumulative DVH Segments Volume');
hold on
plot(DVH(:,3)*1000, (DVH(:,4)*v100_lungs_ml)/100);
hold on
plot(DVH(:,5)*1000, (DVH(:,6)*v100_liver_ml)/100);
hold on
plot(DVH(:,7)*1000, (DVH(:,8)*v100_heart_ml)/100);
hold on
plot(DVH(2:size(DVH,1),9)*1000, (DVH(2:size(DVH,1),10)*v_100_brain_ml)/100);
hold off
legend('Kidneys', 'Lungs', 'Liver', 'Heart', 'Brain');
saveas(all_DVH_cum_ml, 'Cumulative DVH Segments Volume.png');

%Read the DVH file taken from the platform 3D Slicer, using the csvread command,
%from the second line (the first is subtitles). Use of the plot command to
%draw the graphs for each segment: dose in mGy (multiplication needed by 1000)
%vs volume in % and dose in mGy vs volume in ml. In order to obtain the
%volume in ml, it is necessary to multiply by the total volume of the
%segment and divide by 100. 
%Use of the saveas command to save the graphics on the computer. At the end, all plots
%are drawn on the same graph for easy viewing, with legend.

figure;
[kidneys_DVH_plot_dif_frac, max_kidneysf_dDVHy, mean_kidneys_dDVH]=myplot(kidneys_DVH_frac);
ylim([0 max_kidneysf_dDVHy]);
xlabel('Dose (mGy)');
ylabel('Volume Fraction (%)');
title('Differential DVH Kidneys Fraction');
saveas(kidneys_DVH_plot_dif_frac, 'Differential DVH Kidneys Fraction.png');
figure;
[kidneys_DVH_plot_dif_ml, max_kidneysv_dDVH]=myplot(kidneys_DVH_ml);
ylim([0 max_kidneysv_dDVH]);
xlabel('Dose (mGy)');
ylabel('Volume (ml)');
title('Differential DVH Kidneys Volume');
saveas(kidneys_DVH_plot_dif_ml, 'Differential DVH Kidneys Volume.png');
figure;
[lungs_DVH_plot_dif_frac,max_lungsf_dDVHy, mean_lungs_dDVH]=myplot(lungs_DVH_frac);
ylim([0 max_lungsf_dDVHy]);
xlabel('Dose (mGy)');
ylabel('Volume Fraction (%)');
title('Differential DVH Lungs Fraction');
saveas(lungs_DVH_plot_dif_frac, 'Differential DVH Lungs Fraction.png');
figure;
[lungs_DVH_plot_dif_ml, max_lungsv_dDVHy]=myplot(lungs_DVH_ml);
ylim([0 max_lungsv_dDVHy]);
xlabel('Dose (mGy)');
ylabel('Volume (ml)');
title('Differential DVH Lungs Volume');
saveas(lungs_DVH_plot_dif_ml, 'Differential DVH Lungs Volume.png');
figure;
[liver_DVH_plot_dif_frac, max_liverf_dDVHy, mean_liver_dDVH]=myplot(liver_DVH_frac);
ylim([0 max_liverf_dDVHy]);
xlabel('Dose (mGy)');
ylabel('Volume Fraction (%)');
title('Differential DVH Liver Fraction');
saveas(liver_DVH_plot_dif_frac, 'Differential DVH Liver Fraction.png');
figure;
[liver_DVH_plot_dif_ml, max_liverv_dDVHy]=myplot(liver_DVH_ml);
ylim([0 max_liverv_dDVHy]);
xlabel('Dose (mGy)');
ylabel('Volume (ml)');
title('Differential DVH Liver Volume');
saveas(liver_DVH_plot_dif_ml, 'Differential DVH Liver Volume.png');
figure;
[heart_DVH_plot_dif_frac, max_heartf_dDVHy, mean_heart_dDVH]=myplot(heart_DVH_frac);
ylim([0 max_heartf_dDVHy]);
xlabel('Dose (mGy)');
ylabel('Volume Fraction (%)');
title('Differential DVH Heart Fraction');
saveas(heart_DVH_plot_dif_frac, 'Differential DVH Heart Fraction.png');
figure;
[heart_DVH_plot_dif_ml, max_heartv_dDVHy]=myplot(heart_DVH_ml);
ylim([0 max_heartv_dDVHy]);
xlabel('Dose (mGy)');
ylabel('Volume (ml)');
title('Differential DVH Heart Volume');
saveas(heart_DVH_plot_dif_ml, 'Differential DVH Heart Volume.png');
figure;
[brain_DVH_plot_dif_frac, max_brainf_dDVHy, mean_brain_dDVH]=myplot(brain_DVH_frac);
ylim([0 max_brainf_dDVHy]);
xlabel('Dose (mGy)');
ylabel('Volume Fraction (%)');
title('Differential DVH Brain Fraction');
saveas(brain_DVH_plot_dif_frac, 'Differential DVH Brain Fraction.png');
figure;
[brain_DVH_plot_dif_ml, max_brainv_dDVHy]=myplot(brain_DVH_ml);
ylim([0 max_brainv_dDVHy]);
xlabel('Dose (mGy)');
ylabel('Volume (ml)');
title('Differential DVH Brain Volume');
saveas(brain_DVH_plot_dif_ml, 'Differential DVH Brain Volume.png');

mean_dDVH=[mean_liver_dDVH; mean_heart_dDVH; mean_kidneys_dDVH; mean_lungs_dDVH; mean_brain_dDVH];

table_dDVH=table(organ, mean_dDVH);
writetable(table_dDVH, 'dose_organ_DVH.xlsx', 'Sheet', 2);

figure;
all_DVH_dif_frac=myplot(kidneys_DVH_frac);
ylim([0 max_kidneysf_dDVHy]);
xlabel('Dose (mGy)');
ylabel('Volume Fraction (%)');
title('Differential DVH Segments Fraction');
hold on
myplot(lungs_DVH_frac);
hold on
myplot(liver_DVH_frac);
hold on
myplot(heart_DVH_frac);
hold on
myplot(brain_DVH_frac);
hold off
legend('Kidneys', 'Lungs', 'Liver', 'Heart', 'Brain');
saveas(all_DVH_dif_frac, 'Differential DVH Segments Fraction.png');

figure;
all_DVH_dif_ml=myplot(kidneys_DVH_ml);
ylim([0 max_lungsv_dDVHy]);
xlabel('Dose (mGy)');
ylabel('Volume (ml)');
title('Differential DVH Segments Volume');
hold on
myplot(lungs_DVH_ml);
hold on
myplot(liver_DVH_ml);
hold on
myplot(heart_DVH_ml);
hold on
myplot(brain_DVH_ml);
hold off
legend('Kidneys', 'Lungs', 'Liver', 'Heart', 'Brain');
saveas(all_DVH_dif_ml, 'Differential DVH Segments Volume.png');

%Since the Slicer's dose volume histogram is a cumulative histogram, an
%external script from Mathworks is used to convert the cumulative histogram
%into a differential histogram. Subsequently, the steps are similar to
%those described above.

%In terms of absolute volume, the lungs have greater volume. In terms of
%percentage, it is the kidneys.