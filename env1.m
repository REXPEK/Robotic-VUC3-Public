function [env,ctrl] = env1

% pusher
try
    env.robot = loadrobot('universalUR5','DataFormat','column');
catch ME %sometimes the user might not have the Robotics System Toolbox
    disp('An error occurred while loading the robot (this might indicate a missing Robotics System Toolbox, 3D animation will be disabled):');
    disp(ME.message); 
    env.robot = [];
end

% geometry
env.geom.type = 1;
env.geom.a = .2;
env.geom.b = .1;
env.geom.r = 0.01;

env = dynamics(env);

% initial configuration
env.initial.x = .25;
env.initial.y = .1;
env.initial.t = 0;
env.initial.c = 0;

% goal configuration
env.goal.x = -.25;
env.goal.y = .75;

% time
env.time.dt = .1;
env.time.T = 20;

% simulation
env.sim.b = env.geom.b;

env.sim.lamda_xx = 1.1;
env.sim.lamda_xy = 0;
env.sim.lamda_xt = .1;
env.sim.lamda_yy = .8;
env.sim.lamda_yt = .1;
env.sim.lamda_tt = 1/env.model.beta^2*10;

% noise
env.noise.initial.covariance = diag([1e-3 1e-3 1e-2 1e-3]);
env.noise.f.covariance = diag([1e-1 1e-1 1e-1 1e-1]);
env.noise.u.covariance = diag([1e-1 2]);

env.noise.initial.covariance = diag([0 0 0 0]);
env.noise.f.covariance = diag([0 0 0 0]);
env.noise.u.covariance = diag([0 0]);

% type
ctrl.type = 2;

end

