%EUDMODEL(DVH), where DVH is a 2 column matrix corresponding to the cumulative, not
%differential, dose volume histogram. The 1st column corresponds to increasing absolute dose or 
%percentage dose values, and the 2nd column to the corresponding absolute or relative volume value.
%The matrix must have a minimum of two rows, and both columns must be of equal length.
%by Hiram A. Gay, MD
%Revised July 8 2007
    function [diffDVH, max_dDVH_y, mean_dDVH] = myplot(dvh)
   %verifying that the cumulative DVH has at least 2 rows and columns
   [nb,N]=size(dvh);
  if (nb < 2)
      disp('Error: Cumulative dvh must have at least 2 rows.'); return;
  end
  if (N < 2)
      disp('Error: Cumulative dvh must have at least 2 columns.'); return;
  end
  %verifying that the cumulative DVH has no negative numbers in the dose or volume columns
  for i=1:nb
      if (dvh(i,1) < 0)
          message = sprintf('Error: Dose data error. dvh column 1, row %g is negative',i);
          fprintf(message); return;
      end
      if (dvh(i,2) < 0)
          message = sprintf('Error: Volume data error. dvh column 2, row %g is negative',i);
          fprintf(message); return;
      end
  end
  % Converting cumulative DVH to differential DVH, and checking for DVH errors
  for i=2:nb
      dvh(i-1,1)=dvh(i-1,1)+(dvh(i,1)-dvh(i-1,1))/2;
      if (dvh(i,1)-dvh(i-1,1) <= 0)
          message = sprintf('Error: Dose data error. dvh column 1, row %g <= dvh column 1, row %g',i,i-1);
          fprintf(message); return;
      end
      dvh(i-1,2)=(dvh(i-1,2)-dvh(i,2));
      if (dvh(i-1,2) < 0)
          message = sprintf('Error: Volume bin < 0. Verify dvh column 2, rows %g and %g',i-1,i);
          fprintf(message); return;
      end
  end
   diffDVH=plot(dvh(:,1), dvh(:,2));
   
   max_dDVH_y=max(dvh(:,2));
 
  s=0;
  for i=1:nb
      s=s+(dvh(i,1)*dvh(i,2));
  end
    mean_dDVH=s/100;
    
%     min_dDVH=min(dvh(:,1));

    

end