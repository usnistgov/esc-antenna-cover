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
  
    result = py.circlecover.circles.min_cover_greedy(pcenters,contour,distance);

    ccles = result{1};
    centers_x = zeros(1,length(ccles));
    centers_y = zeros(1,length(ccles));
    radius = zeros(1,length(ccles));
    
    % stuff the return values into a matlab array.
    % There's probably a better way to do this than looping.
    for k = 1:length(ccles)
        ccle = ccles{k};
        centers_x(k) = ccle.center{1};
        centers_y(k) = ccle.center(2);
        radius(k) = ccle.get_radius();
    end
end
        
    
        
        
    