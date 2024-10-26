function [U,dY] = cascadedP(X,env,ctrl)
dxd = ctrl.par(1)*(env.goal.x-X(1));
dyd = ctrl.par(2)*(env.goal.y-X(2));
vd = (env.model.beta^2+X(4)^2)/(env.model.beta^2)*sqrt(dxd^2+dyd^2);
td = -atan2(dxd,dyd);
dtd = ctrl.par(3)*(td-X(3));
cd = fsolve(@(c) dtd*(env.model.beta^2+c^2) - c*vd,X(4),optimset('Display','off'));
dcd = ctrl.par(4)*(cd-X(4));
wd = dcd + env.geom.b/2*X(4)/(env.model.beta^2+X(4)^2)*vd;
U = [vd;wd];
dY = zeros(4,1);
end
