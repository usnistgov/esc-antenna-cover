function  [centers_x,centers_y,radius] = min_cover(centers,ic, distance)
% find the min area greedy cover given the possible centers
% and interference contour.
%
% parameters:
%
%
% ic : interference contour Kx2 matrix giving the points defining the
%      interference contour.
% distance: minimum placement separation
%
% returns:
%
% centers_x: vector giving the x coord of centers where circles 
%    should be placed for the min cover.
% centers_y : vector giving the y-coord of centers where circles should
%    be placed for the min cover.
%
                
    pcenters = py.list();
    for k = 1:length(centers)
        cent = centers(k,:);
        pycent = py.list();
        pycent.append(cent(1));
        pycent.append(cent(2));
        pcenters.append(pycent);
    end;
    
    contour = py.list();
    
    
    for m = 1:length(ic)
        point = ic(m,:);
        pypoint = py.list();
        pypoint.append(point(1));
        pypoint.append(point(2));
        contour.append(pypoint);
    end;
  
    result = py.circlecover.circlecover.min_area_cover_greedy(pcenters,contour,distance);

    ccles = result{1};
    centers_x = zeros(1,length(ccles));
    centers_y = zeros(1,length(ccles));
    radius = zeros(1,length(ccles));
    
    % stuff the return values into a matlab array.
    % There's probably a better way to do this than looping.
    for k = 1:length(ccles)
        ccle = ccles{k};
        centers_x(k) = ccle.center{1};
        centers_y(k) = ccle.center{2};
        radius(k) = ccle.get_radius();
    end
end
        

% This software was developed by employees of the National Institute
% of Standards and Technology (NIST), an agency of the Federal
% Government. Pursuant to title 17 United States Code Section 105, works
% of NIST employees are not subject to copyright protection in the United
% States and are considered to be in the public domain. Permission to freely
% use, copy, modify, and distribute this software and its documentation
% without fee is hereby granted, provided that this notice and disclaimer
% of warranty appears in all copies.
%
% THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND,
% EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
% TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY
% IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
% AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION
% WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE
% ERROR FREE. IN NO EVENT SHALL NASA BE LIABLE FOR ANY DAMAGES, INCLUDING,
% BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES,
% ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS
% SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR
% OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY
% OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT
% OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.
%
% Distributions of NIST software should also include copyright and licensing
% statements of any third-party software that are legally bundled with
% the code in compliance with the conditions of those licenses.
%    
        
    
