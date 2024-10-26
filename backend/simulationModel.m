function dX = simulationModel(X,U,env)

b = env.sim.b;

l_xx = env.sim.lamda_xx;
l_xy = env.sim.lamda_xy;
l_xt = env.sim.lamda_xt;
l_yy = env.sim.lamda_yy;
l_yt = env.sim.lamda_yt;
l_tt = env.sim.lamda_tt;

t = X(3);
c = X(4);

vn = U(1);
vt = U(2);

dx = vn*cos(t)*(l_xt*cos(t)*tan(t)*b - 2*l_xt*tan(t)*sin(t)*c - 2*l_xt*cos(t)*c - l_xt*sin(t)*b + 2*tan(t)*l_xx - 2*l_xy)/(2*(2*cos(t)*l_xy*sin(t) - 2*cos(t)*c*l_yt - l_xx*sin(t)^2 + l_yy*sin(t)^2 + 2*c*l_xt*sin(t) - c^2*l_tt - l_yy));
dy = vn*cos(t)*(l_yt*cos(t)*tan(t)*b - 2*l_yt*tan(t)*sin(t)*c - 2*cos(t)*c*l_yt - l_yt*sin(t)*b + 2*tan(t)*l_xy - 2*l_yy)/(2*(2*cos(t)*l_xy*sin(t) - 2*cos(t)*c*l_yt - l_xx*sin(t)^2 + l_yy*sin(t)^2 + 2*c*l_xt*sin(t) - c^2*l_tt - l_yy));
dt = vn*cos(t)*(l_tt*cos(t)*tan(t)*b - 2*l_tt*tan(t)*sin(t)*c - 2*l_tt*cos(t)*c - l_tt*sin(t)*b + 2*tan(t)*l_xt - 2*l_yt)/(2*(2*cos(t)*l_xy*sin(t) - 2*cos(t)*c*l_yt - l_xx*sin(t)^2 + l_yy*sin(t)^2 + 2*c*l_xt*sin(t) - c^2*l_tt - l_yy));
dc = -(4*l_yy*vt - 4*cos(t)^3*c*l_xt*vn + 4*tan(t)*cos(t)^2*l_xx*vn - 8*sin(t)*cos(t)*l_xy*vt - 4*sin(t)*cos(t)*l_yy*vn - 8*sin(t)*c*l_xt*vt - 2*cos(t)*b*l_yt*vn + 8*cos(t)*c*l_yt*vt - 4*tan(t)*sin(t)^2*cos(t)*c*l_yt*vn + 2*tan(t)*sin(t)*cos(t)^2*b*l_yt*vn - 4*tan(t)*sin(t)*cos(t)^2*c*l_xt*vn - 2*tan(t)*sin(t)*cos(t)*b*c*l_tt*vn + 4*sin(t)^2*l_xx*vt - 4*sin(t)^2*l_yy*vt + 4*c^2*l_tt*vt - 4*cos(t)^2*l_xy*vn + 2*tan(t)*cos(t)^3*b*l_xt*vn + tan(t)*cos(t)^2*b^2*l_tt*vn - 2*sin(t)^2*cos(t)*b*l_yt*vn - 2*sin(t)*cos(t)^2*b*l_xt*vn - 4*sin(t)*cos(t)^2*c*l_yt*vn - sin(t)*cos(t)*b^2*l_tt*vn - 2*cos(t)^2*b*c*l_tt*vn + 4*tan(t)*sin(t)*cos(t)*l_xy*vn + 2*tan(t)*cos(t)*b*l_xt*vn)/(4*(2*cos(t)*l_xy*sin(t) - 2*cos(t)*c*l_yt - l_xx*sin(t)^2 + l_yy*sin(t)^2 + 2*c*l_xt*sin(t) - c^2*l_tt - l_yy));

dX = [dx;dy;dt;dc];

end