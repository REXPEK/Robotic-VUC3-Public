function [U,dY] = cascadedPD(X,Y,env,ctrl)
dY1 = 2/ctrl.par(1)*(-Y(1)) + 1/(ctrl.par(1)^2)*(env.goal.x-X(1));
dY2 = 2/ctrl.par(2)*(-Y(2)) + 1/(ctrl.par(2)^2)*(env.goal.y-X(2));
dxd = Y(1);
dyd = Y(2);
vd = (env.model.beta^2+X(4)^2)/(env.model.beta^2)*sqrt(dxd^2+dyd^2);
td = -atan2(dxd,dyd);
dY3 = 2/ctrl.par(3)*(-Y(3)) + 1/(ctrl.par(3)^2)*(td-X(3));
dtd = Y(3);
cd = fsolve(@(c) dtd*(env.model.beta^2+c^2) - c*vd,X(4),optimset('Display','off'));
dY4 = 2/ctrl.par(4)*(-Y(4)) + 1/(ctrl.par(4)^2)*(cd-X(4));
dcd = Y(4);
wd = dcd + env.geom.b/2*X(4)/(env.model.beta^2+X(4)^2)*vd;
U = [vd;wd];
dY = [dY1;dY2;dY3;dY4];
end
