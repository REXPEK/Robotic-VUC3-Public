function env = dynamics(env)
a = env.geom.a;
b = env.geom.b;
if nargin < 2
    env.geom.type = 1;
end
switch env.geom.type
    case 1
        env.model.beta = (log(sqrt(a^2 + b^2) + b)*a^3 - log(a)*a^3 + log(b)*b^3 - log(sqrt(a^2 + b^2) - a)*b^3 + 2*sqrt(a^2 + b^2)*a*b)/(12*a*b);
    case 2
        env.model.beta = sqrt((a^2+b^2)/12);
end
dx = @(c,v,w) 0;
dy = @(c,v,w) env.model.beta^2/(env.model.beta^2 + c^2)*v;
dt = @(c,v,w) c*v/(env.model.beta^2 + c^2);
dc = @(c,v,w) w -  env.geom.b/2*c/( env.model.beta^2 + c^2)*v;
env.model.f = @(x,u) [rot(x(3))*[dx(x(4),u(1),u(2));
    dy(x(4),u(1),u(2))];
    dt(x(4),u(1),u(2));
    dc(x(4),u(1),u(2))];
end