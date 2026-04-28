% igf not filling

function odata = gapfillclosest( data, igf )

if nargin == 1
  igf = [];
end

[nx,ny,nz] = size(data);

%= saving original data
odata = data;

for z = 1:nz

  idata = squeeze( data(:,:,z) );

  % dum vallue so nan values are not searched
  idata(igf) = nanmean(idata(:));

  [indx,indy] = find(isnan(idata)==1 );
  ind = find( indx > 1 & indx < nx ); 
  indx = indx(ind);
  indy = indy(ind);

  ind = find( indy > 1 & indy < ny ); 
  indx = indx(ind);
  indy = indy(ind);

  na  = length( indx );


  for a = 1:na
    xa = indx(a);
    ya = indy(a);
    idata(xa,ya)  = nanmean([ idata(xa+1,ya) idata(xa-1,ya) idata(xa, ya+1) idata(xa, ya-1) ]);
  end

  [indx,indy] = find(isnan(idata)==1 );

  na  = length( indx );

  for a = 1:na
    
    xa = indx(a);
    ya = indy(a);

    if xa == 1 & ya == 1
      idata(xa,ya)  = nanmean([ idata(xa+1,ya) idata(xa, ya+1)]);
    elseif xa == 1 & ya == ny
      idata(xa,ya)  = nanmean([ idata(xa+1,ya) idata(xa, ya-1) ]);
    elseif  xa == nx & ya == 1
      idata(xa,ya)  = nanmean([ idata(xa-1,ya) idata(xa, ya+1) ]);
    elseif  xa == nx & ya == ny
      idata(xa,ya)  = nanmean([ idata(xa-1,ya) idata(xa, ya-1) ]);
    elseif xa == 1
      idata(xa,ya)  = nanmean([ idata(xa+1,ya) idata(xa, ya+1) idata(xa, ya-1) ]);
    elseif xa == nx
      idata(xa,ya)  = nanmean([ idata(xa-1,ya) idata(xa, ya+1) idata(xa, ya-1) ]);
    elseif ya == 1
      idata(xa,ya)  = nanmean([ idata(xa+1,ya) idata(xa-1,ya) idata(xa, ya+1) ]);
    elseif ya == ny
      idata(xa,ya)  = nanmean([ idata(xa+1,ya) idata(xa-1,ya) idata(xa, ya-1) ]);
    end

  end

  ni = length( find( isnan(idata) == 1 ));
  if ni > 0
    disp( [ 'Z=', num2str(z), ': Filling ', num2str(na), ' cells, ', num2str(100*na/nx/ny), '% of total, ', num2str(ni), ' left.'] );
  end

  if nargin == 2
    idata(igf) = nan;
  end
  odata(:,:,z) = idata;



end

   


return