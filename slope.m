% determine slope angle and slope aspect
%source: http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-slope-works.htm
%        http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-aspect-works.htm 

%sample
data=[101 92 85;101 92 85; 101 91 84];
cellsize=1;

dz_dx=((data(1,3)+(2*data(2,3))+data(3,3))-(data(1,1)+(2*data(2,1))+data(3,1)))/8*cellsize;
dz_dy=((data(3,1)+(2*data(3,2))+data(3,3))-(data(1,1)+(2*data(1,2))+data(1,3)))/8*cellsize;

aspect=57.29578*atan2(dz_dy,-dz_dx); % 1 rad = 57.29578
rise_run=sqrt((dz_dx^2)+(dz_dy^2));
slope_degrees=atan(rise_run)*57.29578;
