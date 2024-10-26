function vis3(x,y,theta,a,b,color)
height = .01;
R = rot(theta);
p = diag([a/2,b/2])*[1 1 -1 -1 1;-1 1 1 -1 -1];
pp = R*p + [x;y];
fill3(pp(1,:),pp(2,:),zeros(size(pp(1,:))),color,'edgecolor','none'); %set(h,'facealpha',.2);
fill3(pp(1,:),pp(2,:),height*ones(size(pp(1,:))),color,'edgecolor','none'); %set(h,'facealpha',.2);
z = height*[1 1 0 0 1];
for i = 1:4
    fill3(pp(1,[0 1 1 0 0]+i),pp(2,[0 1 1 0 0]+i),z,color,'edgecolor','none'); %set(h,'facealpha',.2);
end
end