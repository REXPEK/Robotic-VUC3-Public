function vis(x,y,theta,a,b,color)
R = rot(theta);
p = diag([a/2,b/2])*[1 1 -1 -1 1;-1 1 1 -1 -1];
pp = R*p + [x;y];
h = fill(pp(1,:),pp(2,:),color,'edgecolor','k'); set(h,'facealpha',.1);
end