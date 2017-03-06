function  [centers_x,centers_y,indexes,azimuth_angles] = min_isotropic_area_cover(centers, interference_contour, coverage_file, antenna_angle, distance)
% find the min area greedy cover given the possible centers
% and interference contour.
%
% Parameters:
%
% centers: A Mx2 matrix giving the possible locations where sensors can be placed (i.e. the 
%    locations on the shore where antennas can be placed).
% interference contour: Kx2 matrix giving the points defining the
%      interference contour.
% coverage_file: Antenna detection coverage file. This is a set of concentric 
%           bounding regions arranged from smallest to largest.
%           This is specified as a list of one dimensional vectors.
%           The first vector is a list of antenna orientations.
%           The remaining vectors define coverage.
%           Each coverage vector consists of set of coordinates defining the antenna coverage for 
%           the orientation specified in the first vector, 
%           assuming that the antennas have been placed at location (0,0).
%           The coverages can be generated using propagation modeling.
% antenna_angle: 3dB angle for antenna aperture
% distance: minimum placement separation.
%
%
% Returns:
%
% centers_x: vector giving the x coord of centers where circles 
%    should be placed for the min cover.
% centers_y : vector giving the y-coord of centers where circles should
%    be placed for the min cover.
% indexes: A vector that gives the index into the coverage file. For example
%           the index 6 would indicate the 6th index of the coverages specifed 
%           by the coverage_file in the parameter list.
% angles: The azimuth angle of the antenna corresponding to the indexes vector.
%
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
        point = interference_contour(m,:);
        pypoint = py.list();
        pypoint.append(point(1));
        pypoint.append(point(2));
        contour.append(pypoint);
    end;
  
    result = py.antennacover.min_antenna_area_cover_anneal(pcenters, contour, coverage_file, antenna_angle,  min_center_distance)

    rcenters = result{1};
    % 1 added to support MATLAB indexing.
    indexes = result{2} + 1;
    azimuth_angles = result{3};
    
    centers_x = zeros(1,length(rcenters));
    centers_y = zeros(1,length(rcenters));
    
    % stuff the return values into a matlab array.
    % There's probably a better way to do this than looping.
    for k = 1:length(rcenters)
        centers_x(k) = rcenters[k]{1};
        centers_y(k) = rcenters[k]{2};
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
        
    
