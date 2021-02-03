function create3Ddisc(center, radius, borderSize,color,parent)

theta=0:0.01:2*pi;
v=null([0, 0, 1]); %normal
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

plot3(points(1,:),points(2,:),points(3,:),'-','Color',color,'LineWidth',borderSize,'parent',parent);


end