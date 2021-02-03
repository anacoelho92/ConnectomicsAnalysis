function [x,y,z]=sphere2(a,b,c,R,C,transp,parent)
% 
% P1 - center of the sphere
% R - radius of cylinder
% C = colour
 

r=R;

phi=linspace(0,pi,30);
theta=linspace(0,2*pi,40);
[phi,theta]=meshgrid(phi,theta);

x=r*sin(phi).*cos(theta);
y=r*sin(phi).*sin(theta);
z=r*cos(phi); 

if size(C,2)==3
    colour=C;
elseif size(C,2)==1
    switch C
        case 1
            colour='r';
        case 2
            colour='b';
        otherwise
            colour=[rand rand rand];
    end
else
    error('NO SPHERE COLOR IEDNTIFIED.');
end

if nargin==7 && ~isempty(parent)
    sph=surf(x+a, y+b, z+c,'FaceColor',colour,'EdgeColor','none', 'parent', parent);
else
    sph=surf(x+a, y+b, z+c,'FaceColor',colour,'EdgeColor','none');
end

if nargin >= 6 && ~isempty(transp) && transp > 0 && transp < 1
    set(sph,'FaceAlpha', transp)  
end


% shading interp
% camlight
% grid on
 
